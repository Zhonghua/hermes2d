// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "common.h"
#include "limit_order.h"
#include "linsystem.h"
#include "nonlinsystem.h"
#include "weakform.h"
#include "solver.h"
#include "space.h"
#include "precalc.h"
#include "refmap.h"
#include "solution.h"
#include "integrals_h1.h"
#include "views/view.h"
#include "views/vector_view.h"

#include "python_solvers.h"

void NonlinSystem::init_nonlin()
{
  alpha = 1.0;
  res_l2 = res_l1 = res_max = -1.0;

  // Tell LinSystem not to add Dirichlet contributions to the RHS.
  // The reason for this is that in NonlinSystem the Jacobian matrix 
  // is assembled, and the Dirichlet lift is cancelled by the derivative 
  // with respect to the coefficient vector.
  want_dir_contrib = false;
}

// this is needed because of a constructor in RefSystem
NonlinSystem::NonlinSystem() {}

NonlinSystem::NonlinSystem(WeakForm* wf_, 
    Solver* solver_)
{ 
  this->init_lin(wf_, solver_);
  this->init_nonlin();
}

NonlinSystem::NonlinSystem(WeakForm* wf_)
{
  Solver *solver_ = NULL;
  this->init_lin(wf_, solver_);
  this->init_nonlin();
}

NonlinSystem::NonlinSystem(WeakForm* wf_, Solver* solver_, 
    Tuple<Space*> spaces_)
{
  int n = spaces_.size();
  if (n != wf_->neq) 
    error("Number of spaces does not match number of equations in LinSystem::LinSystem().");
  this->init_lin(wf_, solver_);
  this->init_spaces(spaces_);
  this->alloc_vectors();
  this->init_nonlin();
}

NonlinSystem::NonlinSystem(WeakForm* wf_, 
    Tuple<Space*> spaces_)
{
  Solver* solver_ = NULL;
  this->init_lin(wf_, solver_);
  this->init_spaces(spaces_);
  this->alloc_vectors();
  this->init_nonlin();
}

NonlinSystem::NonlinSystem(WeakForm* wf_, Solver* solver_, 
    Space *s_) 
{
  if (wf_->neq != 1) 
    error("Number of spaces does not match number of equations in LinSystem::LinSystem().");
  this->init_lin(wf_, solver_);
  this->init_space(s_);
  this->alloc_vectors();
  this->init_nonlin();
}

NonlinSystem::NonlinSystem(WeakForm* wf_, 
    Space *s_)
{
  Solver *solver_ = NULL;
  this->init_lin(wf_, solver_);
  this->init_space(s_);
  this->alloc_vectors();
  this->init_nonlin();
}

void NonlinSystem::free()
{
  /* FIXME - MEMORY LEAK
  LinSystem::free_matrix();
  LinSystem::free_vectors();
  if (solver) solver->free_data(slv_ctx);

  struct_changed = values_changed = true;
  memset(sp_seq, -1, sizeof(int) * wf->neq);
  wf_seq = -1;
  */
}

void NonlinSystem::assemble(bool rhsonly)
{
  int ndof = this->get_num_dofs();
  if (rhsonly) error("Parameter rhsonly = true has no meaning in NonlinSystem.");

  // assemble J(Y_n) and store in A, assemble F(Y_n) and store in RHS
  LinSystem::assemble();

  // calculate norms of the residual F(Y_n)
  res_l2 = res_l1 = res_max = 0.0;
  for (int i = 0; i < ndof; i++)
  {
    res_l2 += sqr(RHS[i]);
    res_l1 += magn(RHS[i]);
    if (magn(RHS[i]) > res_max) res_max = magn(RHS[i]);
  }
  res_l2 = sqrt(res_l2);

  // multiply RHS by -alpha
  for (int i = 0; i < ndof; i++)
    RHS[i] *= -alpha;
}


bool NonlinSystem::solve(Tuple<Solution*> sln)
{
  // if the number of solutions does not match the number of equations, throw error
  int n = sln.size();
  if (n != this->wf->neq) 
    error("Number of solutions does not match the number of equations in LinSystem::solve().");

  // if no matrix solver defined, throw error
  if (!this->solver) error("No matrix solver defined in NonlinSystem::solve().");

  // if Vec is not initialized, throw error
  if (this->Vec == NULL) error("Vec is NULL in NonlinSystem::solve().");

  // if vector length is not equal to matrix size, throw error
  int ndof = this->get_num_dofs();
  if (ndof != this->A->get_size()) 
    error("Matrix size does not match vector length in NonlinSystem:solve().");

  // The solve() function is almost identical to the original one in LinSystem
  // except that Y_{n+1} = Y_{n} + dY_{n+1}
  TimePeriod cpu_time;

  // solve the system - this is different from LinSystem
  scalar* delta = (scalar*) malloc(ndof * sizeof(scalar));
  memcpy(delta, RHS, sizeof(scalar) * ndof);
  solve_linear_system_scipy_umfpack(this->A, delta);
  report_time("Solved in %g s", cpu_time.tick().last());
  // add the increment dY_{n+1} to the previous solution vector
  for (int i = 0; i < ndof; i++) Vec[i] += delta[i];
  ::free(delta);

  // copy the solution coefficient vectors into Solutions
  cpu_time.tick(H2D_SKIP);
  for (int i = 0; i < n; i++)
  {
    sln[i]->set_fe_solution(spaces[i], pss[i], Vec);
  }
  report_time("Exported solution in %g s", cpu_time.tick().last());
  return true;
}

// single equation case
bool NonlinSystem::solve(Solution* sln)
{
  bool flag;
  flag = this->solve(Tuple<Solution*>(sln));
  return flag;
}

// two equations case
bool NonlinSystem::solve(Solution* sln1, Solution* sln2)
{
  bool flag;
  flag = this->solve(Tuple<Solution*>(sln1, sln2));
  return flag;
}

// three equations case
bool NonlinSystem::solve(Solution* sln1, Solution* sln2, Solution* sln3)
{
  bool flag;
  flag = this->solve(Tuple<Solution*>(sln1, sln2, sln3));
  return flag;
}

// Newton's loop for one equation
bool NonlinSystem::solve_newton(Solution* u_prev, double newton_tol, int newton_max_iter,
                                  Filter* f1, Filter* f2, Filter* f3) {
    int it = 1;
    double res_l2_norm;
    Solution sln_iter;
    Space *space = this->get_space(0);
    do
    {
      info("---- Newton iter %d:", it); it++;
      //printf("ndof = %d\n", space->get_num_dofs());

      // reinitialization of filters (if relevant)
      if (f1 != NULL) f1->reinit();
      if (f2 != NULL) f2->reinit();
      if (f3 != NULL) f3->reinit();

      // assemble the Jacobian matrix and residual vector,
      // solve the system
      this->assemble();
      this->solve(&sln_iter);

      // calculate the l2-norm of residual vector
      res_l2_norm = this->get_residual_l2_norm();
      info("Residual L2 norm: %g", res_l2_norm);

      // save the new solution as "previous" for the
      // next Newton's iteration
      u_prev->copy(&sln_iter);
    }
    while (res_l2_norm > newton_tol && it <= newton_max_iter);
    
    // returning "true" if converged, otherwise returning "false"
    if (it <= newton_max_iter) return true;
    else return false;
}

// Newton's loop for two equations
bool NonlinSystem::solve_newton(Solution* u_prev_1, Solution* u_prev_2, double newton_tol, int newton_max_iter,
                                  Filter* f1, Filter* f2, Filter* f3) {
    int it = 1;
    double res_l2_norm;
    Solution sln_iter_1, sln_iter_2;
    Space *space_1 = this->get_space(0);
    Space *space_2 = this->get_space(1);
    do
    {
      info("---- Newton iter %d:", it); it++;
      //printf("ndof = %d\n", space_1->get_num_dofs() + space_2->get_num_dofs());

      // reinitialization of filters (if relevant)
      if (f1 != NULL) f1->reinit();
      if (f2 != NULL) f2->reinit();
      if (f3 != NULL) f3->reinit();

      // assemble the Jacobian matrix and residual vector,
      // solve the system
      this->assemble();
      this->solve(Tuple<Solution*>(&sln_iter_1, &sln_iter_2));

      // calculate the l2-norm of residual vector
      res_l2_norm = this->get_residual_l2_norm();
      info("Residual L2 norm: %g", res_l2_norm);

      // save the new solutions as "previous" for the
      // next Newton's iteration
      u_prev_1->copy(&sln_iter_1);
      u_prev_2->copy(&sln_iter_2);
    }
    while (res_l2_norm > newton_tol && it <= newton_max_iter);

    // returning "true" if converged, otherwise returning "false"
    if (it <= newton_max_iter) return true;
    else return false;
}

// Newton's loop for three equations
bool NonlinSystem::solve_newton(Solution* u_prev_1, Solution* u_prev_2, Solution* u_prev_3,
                                  double newton_tol, int newton_max_iter,
                                  Filter* f1, Filter* f2, Filter* f3) {
    int it = 1;
    double res_l2_norm;
    Solution sln_iter_1, sln_iter_2, sln_iter_3;
    Space *space_1 = this->get_space(0);
    Space *space_2 = this->get_space(1);
    Space *space_3 = this->get_space(2);
    do
    {
      info("---- Newton iter %d:", it); it++;

      // reinitialization of filters (if relevant)
      if (f1 != NULL) f1->reinit();
      if (f2 != NULL) f2->reinit();
      if (f3 != NULL) f3->reinit();

      // assemble the Jacobian matrix and residual vector,
      // solve the system
      this->assemble();
      this->solve(Tuple<Solution*>(&sln_iter_1, &sln_iter_2, &sln_iter_3));

      // calculate the l2-norm of residual vector
      res_l2_norm = this->get_residual_l2_norm();
      //info("Residual L2 norm: %g", res_l2_norm);
      printf("Residual L2 norm: %g\n", res_l2_norm);



  // Show the solution at the end of time step.
  // Initialize views.
  VectorView vview("velocity [m/s]", 0, 0, 500, 400);
  ScalarView pview("pressure [Pa]", 510, 0, 500, 400);
  //vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  pview.show_mesh(true);
  char title[100];
  sprintf(title, "Velocity, iter %d", it);
      vview.set_title(title);
      vview.show(&sln_iter_1, &sln_iter_2, H2D_EPS_LOW);
      sprintf(title, "Pressure, iter %d", it);
      pview.set_title(title);
      pview.show(&sln_iter_3);
      pview.wait(H2DV_WAIT_KEYPRESS);



      // save the new solutions as "previous" for the
      // next Newton's iteration
      u_prev_1->copy(&sln_iter_1);
      u_prev_2->copy(&sln_iter_2);
      u_prev_3->copy(&sln_iter_3);
    }
    while (res_l2_norm > newton_tol && it <= newton_max_iter);

    // returning "true" if converged, otherwise returning "false"
    if (it <= newton_max_iter) return true;
    else return false;
}
