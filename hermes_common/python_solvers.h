// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_PYTHON_SOLVERS_H
#define __HERMES_COMMON_PYTHON_SOLVERS_H

void solve_linear_system_numpy(CooMatrix *mat, double *res);
void solve_linear_system_scipy_umfpack(CooMatrix *mat, double *res);
void solve_linear_system_scipy_cg(CooMatrix *mat, double *res);
void solve_linear_system_scipy_gmres(CooMatrix *mat, double *res);

#endif
