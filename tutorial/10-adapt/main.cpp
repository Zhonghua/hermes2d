#define H2D_REPORT_INFO
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example shows how to run adaptive hp-FEM, h-FEM and p-FEM with
// basic control parameters. The underlying problem is a planar model
// of an electrostatic micromotor (MEMS). You may want to experiment with 
// various types of adaptivity via the options H2D_P_ISO, H2D_P_ANISO, 
// H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H, H2D_HP_ANISO_P, 
// and H2D_HP_ANISO. See the User Documentation for more details. 
//   Uniform initial polynomial degree of mesh elements can be set using
// the variable P_INIT. Before using adaptivity, you have to define a refinement 
// selector as shown below. The function adapt() takes the selector
// as a parameter, along with THRESHOLD, STRATEGY, and MESH_REGULARITY. 
//   Additional control parameters are available, these will be demonstrated
// in the following tutorial examples. In this example, two types of convergence  
// graphs are created -- error estimate wrt. the number of degrees of freedom 
// (DOF), and error estimate wrt. CPU time. Later we will show how to output 
// the error wrt. exact solution when exact solution is available. 
//   This example also demonstrates how to define different material parameters
// in various parts of the computational domain, and how to measure time. 
//
// PDE: -div[eps_r(x,y) grad phi] = 0
//      eps_r = EPS_1 in Omega_1 (surrounding air)
//      eps_r = EPS_2 in Omega_2 (moving part of the motor)
//
// BC: phi = 0 V on Gamma_1 (left edge and also the rest of the outer boundary
//     phi = VOLTAGE on Gamma_2 (boundary of stator)
//
// The following parameters can be changed:

const int P_INIT = 2;                      // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.2;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_H; // Predefined list of element refinement candidates. Possible values are
                                           // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                           // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1; // Maximum allowed level of hanging nodes:
                                // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                // Note that regular meshes are not supported, this is due to
                                // their notoriously bad performance.
const double ERR_STOP = 1.0;    // Stopping criterion for adaptivity (rel. error tolerance between the
const double CONV_EXP = 1.0;    // Default value is 1.0. This parameter influences the selection of
                                // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;    // Adaptivity process stops when the number of degrees of freedom grows
                                // over this limit. This is to prevent h-adaptivity to go on forever.

// Problem parameters.
const int OMEGA_1 = 1;
const int OMEGA_2 = 2;
const int STATOR_BDY = 2;
const double EPS_1 = 1.0;       // Relative electric permittivity in Omega_1.
const double EPS_2 = 10.0;      // Relative electric permittivity in Omega_2.
const double VOLTAGE = 50.0;    // Voltage on the stator.

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return (ess_bdy_marker == STATOR_BDY) ? VOLTAGE : 0.0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("motor.mesh", &mesh);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(biform1), H2D_SYM, OMEGA_1);
  wf.add_matrix_form(callback(biform2), H2D_SYM, OMEGA_2);

  // Initialize views.
  ScalarView sview("Scalar potential Phi", 0, 0, 400, 600);
  VectorView gview("Gradient of Phi", 410, 0, 400, 600);
  gview.set_min_max_range(0, 1e8);
  OrderView  oview("Mesh", 820, 0, 400, 600);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(SOLVER_UMFPACK, space.get_num_dofs(), mat, rhs, solver);

  // Adaptivity loop:
  Solution sln, ref_sln;
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    info("Solving on fine mesh.");

    // Construct the globally refined reference mesh.
    Mesh ref_mesh;
    ref_mesh.copy(&mesh);
    ref_mesh.refine_all_elements();

    // Setup space for the reference solution.
    Space *ref_space = space.dup(&ref_mesh);
    int order_increase = 1;
    ref_space->copy_orders(&space, order_increase);
 
    // Solve the reference problem.
    solve_linear(ref_space, &wf, &ref_sln, SOLVER_UMFPACK);

    // Project the fine mesh solution on the coarse mesh.
    info("Projecting fine mesh solution on coarse mesh.");
    project_global(&space, &ref_sln, &sln);

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution.
    sview.show(&sln);
    gview.show(&sln, &sln, H2D_EPS_NORMAL, H2D_FN_DX_0, H2D_FN_DY_0);
    oview.show(&space);

    // Time measurement.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error.");
    H1Adapt hp(&space);
    hp.set_solutions(&sln, &ref_sln);
    double err_est = hp.calc_error() * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%", 
         space.get_num_dofs(), ref_space->get_num_dofs(), err_est);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(space.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu.dat");

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) done = true;
    else {
      info("Adapting coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (space.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the fine mesh solution - the final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(&ref_sln);
  gview.show(&ref_sln, &ref_sln, H2D_EPS_HIGH, H2D_FN_DX_0, H2D_FN_DY_0);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

