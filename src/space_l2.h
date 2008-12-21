// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@utep.edu>
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

// $Id: space_l2.h 1086 2008-10-21 09:05:44Z jakub $

#ifndef __HERMES2D_SPACE_L2
#define __HERMES2D_SPACE_L2

#include "space.h"


/// L2Space represents a space of scalar functions with discontinuities along
/// mesh edges.
///
///
class L2Space : public Space
{
public:
  
  L2Space(Mesh* mesh, Shapeset* shapeset);
  virtual ~L2Space();

  virtual Space* dup(Mesh* mesh) const;

  virtual int get_edge_order(Element* e, int edge);

protected:
  
  struct L2Data
  {
    int vdof[4];
    int edof[4];
  };
  
  L2Data* ldata;
  int lsize;
  
  virtual void resize_tables();

  virtual void assign_vertex_dofs();
  virtual void assign_edge_dofs();
  virtual void assign_bubble_dofs();

  virtual void get_vertex_assembly_list(Element* e, int iv, AsmList* al);
  virtual void get_edge_assembly_list_internal(Element* e, int ie, AsmList* al);
  virtual void get_bubble_assembly_list(Element* e, AsmList* al);

  static double** l2_proj_mat;
  static double*  l2_chol_p;
  static int      l2_proj_ref;

  virtual scalar* get_bc_projection(EdgePos* ep, int order);
  
};



#endif
