#include "hermes2d.h"
#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int P_INIT = 1;
double EPS = 10e-14;
// This test make sure the Lobatto shape 
// function is linearly independent.
// This is not a complete test, since in this test we only proof the shape function is linearly independent in int P_INIT = 1.

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  // We load the mesh on a (-1, 1)^2 domain.
  mloader.load("ref_square.mesh", &mesh);            

  // Create an H1 space with default shapeset,
  // natural BC, and linear elements.
  H1Space space(&mesh);
  // The type of element, mesh_mode = 4 means a rectangle element.
  int mesh_mode = 4;
  int n = space.get_num_dofs();
  space.set_uniform_order(P_INIT);

  printf("\n.........................\n");

  int *fn_idx = new int [n];
  int m = 0;

  // Get vertex fns index.
  for (int i = 0; i < mesh_mode; i++, m++)
  {
    fn_idx[m] = space.get_shapeset()->get_vertex_index(i);
    printf("m = %d, get_vertex_index(m) = %d\n",
          m, space.get_shapeset()->get_vertex_index(m));
  }
  // Get the edge fns index.
  int order = P_INIT;
  for (int edge_order = 2; edge_order <= order; edge_order++)
  {
    for (int j = 0; j < mesh_mode; j++, m++)
    {
      fn_idx[m] = space.get_shapeset()->get_edge_index(j, 0, edge_order);
      printf("m = %d, get_edge_index(m) = %d\n", m, fn_idx[m]);
    }
  }
  // Get the bubble fns index.
  // Note that : This could not get the number of bubbles by using
  //             space.get_shapeset()->get_num_bubbles(order) or maybe 
  //             I don't know how to get the number of bubbles in a element.
  // If the order of the element greater than 2, then there 
  // should be some bubble functions, and I think we could get 
  // the number of bubbles by using "get_num_bubbles(order)" here. 
  // But it is always return to 0, it is not correct.
  int *bubble_idx = space.get_shapeset()->get_bubble_indices(order);
  for (int i = 0; space.get_shapeset()->get_num_bubbles(order); i++, m++ )
  {
    fn_idx[m] = bubble_idx[i];
    printf("m = %d, get_bubble_index(m) = %d\n", m, fn_idx[m]);
  }
  printf("assembling matrix ...\n");

  // The following code was assume that the order of this element is 1.
  // According to the property of linearly independent, 
  // if a part of the Vector Group is linearly independent, 
  // then the whole Vector Group is linearly independent.
  // Here we proof the four vertex's shape function  is linearly independent. 
  CooMatrix mat(4);
  Vector* rhs = new AVector(4);
  CommonSolverSciPyUmfpack solver;

  printf("Get the four times four matrix\n"); 
  for (int j = 0; j < 4; j++)
  {
    mat.add(0, j, space.get_shapeset()->get_fn_value(fn_idx[j], -1, -1, 0));
    printf("get fn value = %f\n", 
           space.get_shapeset()->get_fn_value(fn_idx[j],-1,-1,0));
  }

  for (int j = 0; j < 4; j++)
  {
    mat.add(1, j, space.get_shapeset()->get_fn_value(fn_idx[j], 1, -1, 0));
    printf("get fn value = %f\n", 
           space.get_shapeset()->get_fn_value(fn_idx[j],1,-1,0));
  }

  for (int j = 0; j < 4; j++)
  {
    mat.add(2, j, space.get_shapeset()->get_fn_value(fn_idx[j], 1, 1, 0));
    printf("get fn value = %f\n", 
           space.get_shapeset()->get_fn_value(fn_idx[j],1,1,0));
  }

  for (int j = 0; j < 4; j++)
  {
    mat.add(3, j, space.get_shapeset()->get_fn_value(fn_idx[j], -1, 1, 0));
    printf("get fn value = %f\n", 
           space.get_shapeset()->get_fn_value(fn_idx[j],-1,1,0));
  }

  printf("Add rhs\n");
  for (int i = 0; i < 4; i++)
  {
    printf("i = %d\n", i);
    rhs->add(i, 0.0);
  }

  solver.solve(&mat, rhs);

  for (int i = 0; i < 4; i++)
  {
    // Get the value of the matrix solution by calling Vertex::get().
    if (rhs->get(i) >= EPS)
    {
      printf("Shape functions are not linearly independent\n");
      return ERROR_FAILURE; 
    }
  }
  printf("Success!\n");
 
  delete rhs;
  delete [] fn_idx;
  return ERROR_SUCCESS;
}

