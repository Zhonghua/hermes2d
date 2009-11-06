#include "hermes2d.h"
#include "solver_umfpack.h"

//  This is a simple non-elliptic problem with known exact solution where one
//  can see the advantages of anisotropic refinements. During each computation, 
//  approximate convergence curves are saved into the files "conv_dof.gp" and 
//  "conv_cpu.gp". As in other adaptivity examples, you can compare hp-adaptivity 
//  (ADAPT_TYPE = 0) with h-adaptivity (ADAPT_TYPE = 1) and p-adaptivity (ADAPT_TYPE = 2).
//  You can turn off and on anisotropic element refinements via the ISO_ONLY
//  parameter. 
//
//  PDE: -div(A(x,y) grad u) = 0
//  where a(x,y) = R in the first and third quadrant
//               = 1 in the second and fourth quadrant 
//
//  Exact solution: u(x,y) = cos(M_PI*y/2)    for x < 0
//                  u(x,y) = cos(M_PI*y/2) + pow(x, alpha)   for x > 0   where alpha > 0
//
//  Domain: square, see the file singpert.mesh
//
//  BC:  Homogeneous Dirichlet
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 1;       // number of initial mesh refinements (the original mesh is just one element)
const int INIT_REF_NUM_BDY = 0;   // number of initial mesh refinements towards the boundary
const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;     // This is a quantitative parameter of the adapt(...) function and
                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;           // Adaptive strategy:
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
const double ERR_STOP = 0.6;      // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;     // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// problem constants
const double R = 161.4476387975881;      // Equation parameter.
const double TAU = 0.1;                  // Equation parameter.
const double RHO = M_PI/4.;              // Equation parameter
const double SIGMA = -14.92256510455152; // Equation parameter

// exact solution
static double fn(double x, double y)
{
  double theta = atan2(y,x);
  if (theta < 0) theta = theta + 2*M_PI;
  double r = sqrt(x*x + y*y);

  double mu;
  if (theta <= M_PI/2.) {
    mu = cos((M_PI/2. - SIGMA)*TAU)*cos((theta-M_PI/2+RHO)*TAU);
  }
  else {
    if (theta <= M_PI) {
      mu = cos(RHO*TAU)*cos((theta - M_PI + SIGMA)*TAU);
    }
    else {
      if (theta <= 3*M_PI/2.) {
        mu = cos(SIGMA*TAU)*cos((theta - M_PI + RHO)*TAU);
      }
      else {
        mu = cos((M_PI/2. - RHO)*TAU)*cos((theta-3*M_PI/2. - SIGMA)*TAU);
      }
    }
  }

  return pow(r, TAU) * mu;
}

static double fndd(double x, double y, double& dx, double& dy)
{
  double theta = atan2(y,x);
  if (theta < 0) theta = theta + 2*M_PI;
  double r = sqrt(x*x + y*y);
  // x-derivative
  if (theta <= M_PI/2.) {
    //mu = cos((M_PI/2. - SIGMA)*TAU) * cos((theta - M_PI/2. + RHO)*TAU);
    dx = TAU*x*pow(r, (2.*(-1 + TAU/2.))) * 
    cos((M_PI/2. - SIGMA)*TAU) *
    cos(TAU*(-M_PI/2. + RHO + theta)) +
    (TAU*y*pow(r, TAU)*cos((M_PI/2. - SIGMA)*TAU) * 
    sin(TAU*(-M_PI/2. + RHO + theta))/(r*r));
  }
  else {
    if (theta <= M_PI) {
      //mu = cos(RHO*TAU)*cos((theta - M_PI + SIGMA)*TAU);
      dx = TAU*x * pow(r, (2.*(-1 + TAU/2.))) * cos(RHO*TAU) * 
      cos(TAU*(-M_PI + SIGMA + theta)) +  
      (TAU*y * pow(r, TAU) * cos(RHO*TAU) * 
      sin(TAU*(-M_PI + SIGMA + theta))/(r*r));
    }
    else {
      if (theta <= 3*M_PI/2.) {
        //mu = cos(SIGMA*TAU) * cos((theta - M_PI + RHO)*TAU);
        dx = TAU*x * pow(r, (2.*(-1 + TAU/2.))) * cos(SIGMA*TAU) * 
        cos(TAU*(-M_PI + RHO + theta)) +  
	(TAU*y * pow(r, TAU) * cos(SIGMA*TAU) * 
	 sin(TAU*(-M_PI + RHO + theta))/(r*r));
      }
      else {
        //mu = cos((M_PI/2. - RHO)*TAU) * cos((theta - 3*M_PI/2. - SIGMA)*TAU);
        dx = TAU*x* pow(r, (2*(-1 + TAU/2.))) * 
        cos((M_PI/2. - RHO)*TAU) * 
	cos(TAU*(-3*M_PI/2. - SIGMA + theta)) +  
	(TAU*y*pow(r, TAU) * cos((M_PI/2. - RHO)*TAU) * 
	sin(TAU*(-3.*M_PI/2. - SIGMA + theta))/(r*r));
      }
    }
  }
  // y-derivative
  if (theta <= M_PI/2.) {
    //mu = cos((M_PI/2. - SIGMA)*TAU) * cos((theta - M_PI/2. + RHO)*TAU);
    dy = TAU*y * pow(r, (2*(-1 + TAU/2.))) * 
      cos((M_PI/2. - SIGMA)*TAU) * 
      cos(TAU*(-M_PI/2. + RHO + theta)) -  
      (TAU * pow(r, TAU) * cos((M_PI/2. - SIGMA)*TAU) * 
       sin(TAU*(-M_PI/2. + RHO + theta))*x/(r*r));
  }
  else {
    if (theta <= M_PI) {
      //mu = cos(RHO*TAU)*cos((theta - M_PI + SIGMA)*TAU);
      dy = TAU*y* pow(r, (2*(-1 + TAU/2.))) * cos(RHO*TAU) * 
	cos(TAU*(-M_PI + SIGMA + theta)) -  
        (TAU * pow(r, TAU) * cos(RHO*TAU) * 
	 sin(TAU*(-M_PI + SIGMA + theta))*x/(r*r));
    }
    else {
      if (theta <= 3*M_PI/2.) {
        //mu = cos(SIGMA*TAU) * cos((theta - M_PI + RHO)*TAU);
        dy = TAU*y * pow(r, (2*(-1 + TAU/2.))) * cos(SIGMA*TAU) * 
	  cos(TAU*(-M_PI + RHO + theta)) -  
	  (TAU * pow(r, TAU) * cos(SIGMA*TAU) * 
	   sin(TAU*(-M_PI + RHO + theta))*x/(r*r));
      }
      else {
        //mu = cos((M_PI/2 - RHO)*TAU) * cos((theta-3*M_PI/2 - SIGMA)*TAU);
        dy = TAU*y * pow(r, (2*(-1 + TAU/2.))) * 
        cos((M_PI/2. - RHO)*TAU) * 
	  cos(TAU*(-3*M_PI/2. - SIGMA + theta)) -
	  (TAU * pow(r, TAU) * cos((M_PI/2. - RHO)*TAU) *
	   sin(TAU*((-3*M_PI)/2. - SIGMA + theta))*x/(r*r));
      }
    }
  }
  return fn(x,y);
}

// boundary condition types
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// function values for Dirichlet boundary conditions
scalar bc_values(int marker, double x, double y)
{
  return fn(x, y);
}

template<typename Real>
Real rhs(Real x, Real y)
{
  return 0;
} 

template<typename Real, typename Scalar>
Scalar bilinear_form_I_III(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return R*int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_II_IV(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return 1*int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  //return fn() * int_v<Real, Scalar>(n, wt, v);
  return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  mesh.load("square_quad.mesh");

  // initial mesh refinement
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_all_elements();

  // initialize the shapeset and the cache
  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate basis functions
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form_I_III), SYM, 0);
  wf.add_biform(0, 0, callback(bilinear_form_II_IV), SYM, 1);
  wf.add_liform(0, callback(linear_form));

  // visualize solution and mesh
  ScalarView sview("Coarse solution", 0, 100, 798, 700);
  OrderView  oview("Polynomial orders", 800, 100, 798, 700);

  // matrix solver
  UmfpackSolver solver;

  // convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_captions("Error Convergence for the Singular Line Problem", "Degrees of Freedom", "Error Estimate [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");
  graph.set_log_y();

  // convergence graph wrt. CPU time
  GnuplotGraph graph_cpu;
  graph_cpu.set_captions("Error Convergence for the Singular Line Problem", "CPU Time", "Error Estimate [%]");
  graph_cpu.add_row("exact error", "k", "-", "o");
  graph_cpu.add_row("error estimate", "k", "--");
  graph_cpu.set_log_x();
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

    // solve the coarse mesh problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln_coarse);

    // time measurement
    cpu += end_time();

    // calculate error wrt. exact solution
    ExactSolution exact(&mesh, fndd);
    double error = h1_error(&sln_coarse, &exact) * 100;
    info("\nExact solution error: %g%%", error);

    // view the solution
    sview.show(&sln_coarse);
    oview.show(&space);

    // time measurement
    begin_time();

    // solve the fine mesh problem
    begin_time();
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &sln_fine);

    // calculate error estimate wrt. fine mesh solution
    H1OrthoHP hp(1, &space);
    double err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;
    info("Estimate of error: %g%%", err_est);

    // add entry to DOF convergence graph
    //graph.add_values(0, pow(space.get_num_dofs(), 1./3), error);
    //graph.add_values(1, pow(space.get_num_dofs(), 1./3), err_est);
    graph.add_values(0, space.get_num_dofs(), error);
    graph.add_values(1, space.get_num_dofs(), err_est);
    graph.save("conv_dof.gp");

    // add entry to CPU convergence graph
    graph_cpu.add_values(0, cpu, error);
    graph_cpu.add_values(1, cpu, err_est);
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

  // show the fine solution - this is the final result
  sview.set_title("Final solution");
  sview.show(&sln_fine);

  // wait for keyboard or mouse input
  View::wait("Waiting for keyboard or mouse input.");
  return 0;
}