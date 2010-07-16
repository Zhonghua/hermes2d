#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//
//  This example shows that the worst thing you can ever do is to approximate
//  smooth parts of solutions with uniform low-order meshes. The exact solution
//  to this Poisson problem is u(x,y) = sin(x)*sin(y), defined in the square
//  (0, pi)x(0, pi).
//
//  PDE: -Laplace u = f.
//
//  Known exact solution, see functions fn() and fndd().
//
//  Domain: square domain (0, pi)x(0, pi), mesh file square_quad.mesh.
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

int P_INIT = 1;                          // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1e-4;            // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.

// Exact solution.
static double fn(double x, double y)
{
  return sin(x)*sin(y);
}

static double fndd(double x, double y, double& dx, double& dy)
{
  dx = cos(x)*sin(y);
  dy = sin(x)*cos(y);
  return fn(x, y);
}

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

// Weak forms.
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real>
Real rhs(Real x, Real y)
{
  return 2*sin(x)*sin(y);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh);

  // Avoid zero ndof situation.
  if (P_INIT == 1) {
    if (is_hp(CAND_LIST)) P_INIT++;
    else mesh.refine_element(0, 0);
  }

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  if (is_p_aniso(CAND_LIST))
    space.set_element_order(0, H2D_MAKE_QUAD_ORDER(P_INIT, P_INIT));

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form), H2D_SYM);
  wf.add_vector_form(callback(linear_form));

  // Initialize views.
  ScalarView sview("Coarse mesh solution", 0, 0, 500, 400);
  OrderView  oview("Coarse mesh", 510, 0, 500, 400);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_dof_exact, graph_cpu_est, graph_cpu_exact;

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
    info("Solving on reference mesh.");

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

    // Project the reference mesh solution on the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    project_global(&space, &ref_sln, &sln);

    // Time measurement.
    cpu_time.tick();

    // Calculate error wrt. exact solution.
    info("Calculating error (exact).");
    ExactSolution exact(&mesh, fndd);
    double err_exact = h1_error(&sln, &exact) * 100;

    // View the solution and mesh.
    sview.show(&sln);
    oview.show(&space);

    // Skip exact error calculation and visualization time. 
    cpu_time.tick(HERMES_SKIP);

    // Calculate error estimate wrt. reference solution.
    info("Calculating error (est).");
    H1Adapt hp(&space);
    hp.set_solutions(&sln, &ref_sln);
    double err_est = hp.calc_error() * 100;

    // Report results.
    info("ndof: %d, ref_ndof: %d, err_est: %g%%, err_exact: %g%%", 
         space.get_num_dofs(), ref_space->get_num_dofs(), err_est, err_exact);

    // Add entries to DOF convergence graphs.
    graph_dof_exact.add_values(space.get_num_dofs(), err_exact);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_dof_est.add_values(space.get_num_dofs(), err_est);
    graph_dof_est.save("conv_dof_est.dat");

    // Add entries to CPU convergence graphs.
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact);
    graph_cpu_exact.save("conv_cpu_exact.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) done = true;
    else {
      info("Adapting the coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      if (space.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}


