#include "hermes2d.h"
#include "solver_umfpack.h"

// This example solves adaptively the electric field in a simplified microwave oven.
// The waves are generated using a harmonic surface current on the right-most edge.
// (Such small cavity is present in every microwave oven). There is a circular
// load located in the middle of the main cavity, defined through a different
// permittivity -- see function in_load(...). One can either use a mesh that is
// aligned to the load via curvilinear elements (ALIGN_MESH = true), or an unaligned
// mesh (ALIGN_MESH = false). Convergence graphs are saved both wrt. the dof number
// and cpu time.
//
// PDE: time-harmonic Maxwell's equations;
//      there is circular load in the middle of the large cavity, whose permittivity
//      is different from the rest of the domain
//
// Domain: square cavity with another small square cavity attached from outside
//         on the right
//
// Meshes: you can either use "oven_load_circle.mesh" containing curved elements
//         aligned with the circular load, or "oven_load_square.mesh" which is not
//         aligned.
//
// BC: perfect conductor on the boundary except for the right-most edge of the small
//     cavity, where a harmonic surface current is prescribed
//
// The following parameters can be changed:

const int P_INIT = 2;             // Initial polynomial degree of all mesh elements.
const bool ALIGN_MESH = true;     // if ALIGN_MESH == true, curvilinear elements aligned with the
                                  // circular load are used, otherwise one uses a non-aligned mesh.
const double THRESHOLD = 0.3;     // This is a quantitative parameter of the adapt(...) function and
                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;           // Adaptive strategy:
                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                  //   error is processed. If more elements have similar errors, refine
                                  //   all to keep the mesh symmetric.
                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                  //   than THRESHOLD times maximum element error.
                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                  //   than THRESHOLD.
                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int ADAPT_TYPE = 0;         // Type of automatic adaptivity:
                                  // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                  // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                  // ADAPT_TYPE = 2 ... adaptive p-FEM.
const bool ISO_ONLY = false;      // Isotropic refinement flag (concerns quadrilateral elements only).
                                  // ISO_ONLY = false ... anisotropic refinement of quad elements
                                  // is allowed (default),
                                  // ISO_ONLY = true ... only isotropic refinements of quad elements
                                  // are allowed.
const int MESH_REGULARITY = -1;   // Maximum allowed level of hanging nodes:
                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                  // Note that regular meshes are not supported, this is due to
                                  // their notoriously bad performance.
const double ERR_STOP = 2.0;      // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 40000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.


// problem constants
const double e_0   = 8.8541878176 * 1e-12;
const double mu_0   = 1.256 * 1e-6;
const double e_r = 1.0;
const double mu_r = 1.0;
const double rho = 3820.0;
const double Cp = 7.531000;
const double freq = 1.0*2450000000.0;
const double omega = 2 * M_PI * freq;
const double c = 1 / sqrt(e_0 * mu_0);
const double kappa  = 2 * M_PI * freq * sqrt(e_0 * mu_0);
const double J = 0.0000033333;

// boundary conditions
int e_bc_types(int marker)
{
  if (marker == 2) return BC_ESSENTIAL; // perfect conductor
  else return BC_NATURAL; // impedance
}

// defining load geometry
bool in_load(double x, double y)
{
  double cx = -0.152994121;
  double cy =  0.030598824;
  double r = 0.043273273;
  if (sqr(cx - x) + sqr(cy - y) < sqr(r)) return true;
  else return false;
}

// gamma as a function of x, y
double gam(int marker, double x, double y)
{
  if (ALIGN_MESH && marker == 1) return 0.03;
  if (!ALIGN_MESH && in_load(x,y)) {
    double cx = -0.152994121;  double cy =  0.030598824;
    double r = sqrt(sqr(cx - x) + sqr(cy - y));
    return (0.03 + 1)/2.0 - (0.03 - 1) * atan(10.0*(r -  0.043273273)) / M_PI;
  }
  return 0.0;
}
double gam(int marker, Ord x, Ord y)
{  return 0.0; }


// relative permittivity as a function of x, y
double er(int marker, double x, double y)
{
  if (ALIGN_MESH && marker == 1) return 7.5;
  if (!ALIGN_MESH && in_load(x,y)) {
    double cx = -0.152994121;  double cy =  0.030598824;
    double r = sqrt(sqr(cx - x) + sqr(cy - y));
    return (7.5 + 1)/2.0 - (7.5 - 1) * atan(10.0*(r -  0.043273273)) / M_PI;
  }
  return 1.0;
}
double er(int marker, Ord x, Ord y)
{  return 1.0; }


// bilinear form
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  complex ikappa = complex(0.0, kappa);
  return 1.0/mu_r * int_curl_e_curl_f<Real, Scalar>(n, wt, u, v) -
         ikappa * sqrt(mu_0 / e_0) * int_F_e_f<Real, Scalar>(n, wt, gam, u, v, e) -
         sqr(kappa) * int_F_e_f<Real, Scalar>(n, wt, er, u, v, e);
}

// surface linear form
template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  complex ii = complex(0.0, 1.0);
  return ii * omega * J * int_v1<Real, Scalar>(n, wt, v); // just second component of v, since J = (0, J)
}

// error calculation
template<typename Real, typename Scalar>
Scalar hcurl_form_kappa(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_curl_e_curl_f<Scalar, Scalar>(n, wt, u, v) + sqr(kappa) * int_e_f<Scalar, Scalar>(n, wt, u, v);
}


int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  if (ALIGN_MESH) mesh.load("oven_load_circle.mesh");
  else  mesh.load("oven_load_square.mesh");

  // initialize the shapeset and the cache
  HcurlShapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  HcurlSpace space(&mesh, &shapeset);
  space.set_bc_types(e_bc_types);
  space.set_uniform_order(P_INIT);

  // enumerate basis functions
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form));
  wf.add_liform_surf(0, callback(linear_form_surf));

  // visualize solution and mesh
  VectorView eview("Electric field",0,0,800, 590);
  OrderView ord("Order", 800, 0, 700, 590);

  // matrix solver
  UmfpackSolver solver;

  // convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_captions("Error Convergence for the Waveguide Problem", "Degrees of Freedom", "Error Estimate [%]");
  graph.add_row("error estimate", "-", "o");
  graph.set_log_y();

  // convergence graph wrt. CPU time
  GnuplotGraph graph_cpu;
  graph_cpu.set_captions("Error Convergence for the Waveguide Problem", "CPU Time", "Error Estimate [%]");
  graph_cpu.add_row("error estimate", "-", "o");
  graph_cpu.set_log_y();

  // adaptivity loop
  int it = 1, ndofs;
  bool done = false;
  double cpu = 0.0;
  Solution sln_coarse, sln_fine;
  do
  {
    info("\n---- Adaptivity step %d ---------------------------------------------\n", it++);

    // time measurement
    begin_time();

    // coarse problem
    LinSystem sys(&wf, &solver);
    sys.set_spaces(1, &space);
    sys.set_pss(1, &pss);
    sys.assemble();
    sys.solve(1, &sln_coarse);

    // time measurement
    cpu += end_time();

    // show real part of the solution
    AbsFilter abs(&sln_coarse);
    eview.set_min_max_range(0, 4e3);
    eview.show(&abs);
    ord.show(&space);

    // time measurement
    begin_time();

    // solve the fine mesh problem
    RefSystem ref(&sys);
    ref.assemble();
    ref.solve(1, &sln_fine);

    // calculate error estimate wrt. fine mesh solution
    HcurlOrthoHP hp(1, &space);
    hp.set_biform(0, 0, callback(hcurl_form_kappa));
    double err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;
    info("Hcurl error estimate: %g%%", hcurl_error(&sln_coarse, &sln_fine) * 100);

    // add entry to DOF convergence graph
    graph.add_values(0, space.get_num_dofs(), err_est);
    graph.save("conv_dof.gp");

    // add entry to CPU convergence graph
    graph_cpu.add_values(0, cpu, err_est);
    graph_cpu.save("conv_cpu.gp");

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY);
      ndofs = space.assign_dofs();
      if (ndofs >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu += end_time();
  }
  while (done == false);
  verbose("Total running time: %g sec", cpu);

  // wait for keyboard or mouse input
  printf("Waiting for keyboard or mouse input.\n");
  View::wait();
  return 0;
}

