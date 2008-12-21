// utilita pro generovani edge shape funkci pro quady v prostoru Hcurl
// edge fce jsou gradienty skalarnich lobato edge fci
// bubble funkce jsou legendrovy fce z sede knihy

#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#include <minmax.h>
#endif

void whitney(int i)
{
    printf
    (
      "static double gradleg_quad_p%d_e1_a_0(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * l0(y);\n"
      "}\n\n"
       "static double gradleg_quad_p%d_e1_a_1(double x, double y)\n"
      "{\n"
      "  return -(Legendre%d(x) * l0(y));\n"
      "}\n\n"
     "static double gradleg_quad_p%d_e1_b(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"

      "static double gradleg_quad_p%d_e1_ax_0(double x, double y)\n"
      "{\n"
      "  return Legendre%dx(x) * l0(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e1_ax_1(double x, double y)\n"
      "{\n"
      "  return -(Legendre%dx(x) * l0(y));\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e1_ay_0(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * dl0(y);\n"
      "}\n\n"
       "static double gradleg_quad_p%d_e1_ay_1(double x, double y)\n"
      "{\n"
      "  return -(Legendre%d(x) * dl0(y));\n"
      "}\n\n"
     "static double gradleg_quad_p%d_e1_bx(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e1_by(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n",
     i,i,i, i,i,   i,i, i,i, i,i,i,i,  i,i
    );


    printf
    (
      "static double gradleg_quad_p%d_e2_a(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e2_b_0(double x, double y)\n"
      "{\n"
      "  return l1(x) * Legendre%d(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e2_b_1(double x, double y)\n"
      "{\n"
      "  return -(l1(x) * Legendre%d(y));\n"
      "}\n\n"

      "static double gradleg_quad_p%d_e2_ax(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e2_ay(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e2_bx_0(double x, double y)\n"
      "{\n"
      "  return dl1(x) * Legendre%d(y);\n"
      "}\n\n"
       "static double gradleg_quad_p%d_e2_bx_1(double x, double y)\n"
      "{\n"
      "  return -(dl1(x) * Legendre%d(y));\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e2_by_0(double x, double y)\n"
      "{\n"
      "  return l1(x) * Legendre%dx(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e2_by_1(double x, double y)\n"
      "{\n"
      "  return -(l1(x) * Legendre%dx(y));\n"
      "}\n\n",
    i,i,i,i,i,    i,i, i,i,i,i,i,i,i,i
    );    


    printf
    (
      "static double gradleg_quad_p%d_e3_a_1(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * l1(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e3_a_0(double x, double y)\n"
      "{\n"
      "  return -(Legendre%d(x) * l1(y));\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e3_b(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"

      "static double gradleg_quad_p%d_e3_ax_1(double x, double y)\n"
      "{\n"
      "  return Legendre%dx(x) * l1(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e3_ax_0(double x, double y)\n"
      "{\n"
      "  return -(Legendre%dx(x) * l1(y));\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e3_ay_1(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * dl1(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e3_ay_0(double x, double y)\n"
      "{\n"
      "  return -(Legendre%d(x) * dl1(y));\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e3_bx(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e3_by(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n",
     i,i,i, i,i,i, i,i,i,  i,i, i,i,   i,i
    );

  printf
  (
      "static double gradleg_quad_p%d_e4_a(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e4_b_1(double x, double y)\n"
      "{\n"
      "  return l0(x) * Legendre%d(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e4_b_0(double x, double y)\n"
      "{\n"
      "  return -(l0(x) * Legendre%d(y));\n"
      "}\n\n"

      "static double gradleg_quad_p%d_e4_ax(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e4_ay(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e4_bx_1(double x, double y)\n"
      "{\n"
      "  return dl0(x) * Legendre%d(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e4_bx_0(double x, double y)\n"
      "{\n"
      "  return -(dl0(x) * Legendre%d(y));\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e4_by_1(double x, double y)\n"
      "{\n"
      "  return l0(x) * Legendre%dx(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%d_e4_by_0(double x, double y)\n"
      "{\n"
      "  return -(l0(x) * Legendre%dx(y));\n"
      "}\n\n",   i,i,i,i,i ,   i,i, i,i,i,i,i,i,i,i
  );
}
void edge_fn(int i, int j)
{
 char c1, c2, c3, c4;
 
 if ((j<2) && (i%2))
 {
  if (j == 1) { c1 = '-'; c2 = ' '; c3 = '-'; c4 = ' ';}
  else {c2 = '-'; c1 = ' '; c3 = ' '; c4 = '-';}
  printf
  (
    "static double gradleg_quad_l%d_l%d_a_0(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_a_1(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );

  printf
  (
    "static double gradleg_quad_l%d_l%d_b_0(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c3, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_b_1(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c4, i, j
  );

  printf
  (
    "static double gradleg_quad_l%d_l%d_ax_0(double x, double y)\n"
    "{\n"
    " return %cd2l%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_ax_1(double x, double y)\n"
    "{\n"
    " return %cd2l%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );
 
  printf
  (
    "static double gradleg_quad_l%d_l%d_bx_0(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c3, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_bx_1(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c4, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_ay_0(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_ay_1(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );

  printf
  (
    "static double gradleg_quad_l%d_l%d_by_0(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * d2l%d(y);\n"
    "}\n\n",
    i, j, c3, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_by_1(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * d2l%d(y);\n"
    "}\n\n",
    i, j, c4, i, j
  );
 }
 else if ((i<2) && (j%2))
 {
  
  if (i == 0) { c1 = '-'; c2 = ' '; c3 = '-'; c4 = ' ';}
  else {c2 = '-'; c1 = ' '; c3 = ' '; c4 = '-';}
  printf
  (
    "static double gradleg_quad_l%d_l%d_a_0(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c3, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_a_1(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c4, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_b_0(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_b_1(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );
 
  printf
  (
    "static double gradleg_quad_l%d_l%d_ax_0(double x, double y)\n"
    "{\n"
    " return %cd2l%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c3, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_ax_1(double x, double y)\n"
    "{\n"
    " return %cd2l%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, c4, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_bx_0(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_bx_1(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_ay_0(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c3, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_ay_1(double x, double y)\n"
    "{\n"
    " return %cdl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, c4, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_by_0(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * d2l%d(y);\n"
    "}\n\n",
    i, j, c1, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_by_1(double x, double y)\n"
    "{\n"
    " return %cl%d(x) * d2l%d(y);\n"
    "}\n\n",
    i, j, c2, i, j
  );

 }
 else
 {
  printf
  (
    "static double gradleg_quad_l%d_l%d_a(double x, double y)\n"
    "{\n"
    " return dl%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_b(double x, double y)\n"
    "{\n"
    " return l%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_ax(double x, double y)\n"
    "{\n"
    " return d2l%d(x) * l%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_bx(double x, double y)\n"
    "{\n"
    " return dl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_ay(double x, double y)\n"
    "{\n"
    " return dl%d(x) * dl%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
  printf
  (
    "static double gradleg_quad_l%d_l%d_by(double x, double y)\n"
    "{\n"
    " return l%d(x) * d2l%d(y);\n"
    "}\n\n",
    i, j, i, j
  );
 }
}

void bubble(int p1, int p2)
{
  printf("/* BUBBLE */\n\n");
  printf("/* BUBBLE ( 1 , 0 ) */\n\n");

  for (int i = 0; i <= p1; i++)
  {
    for (int j = 2; j <= p2 + 1; j++)
    {
     printf(
      "static double gradleg_quad_p%dp%d_b1_a(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * l%d(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%dp%d_b1_b(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n",
      i,j,i,j,  i,j
     );
     printf(
      "static double gradleg_quad_p%dp%d_b1_ax(double x, double y)\n"
      "{\n"
      "  return Legendre%dx(x) * l%d(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%dp%d_b1_ay(double x, double y)\n"
      "{\n"
      "  return Legendre%d(x) * dl%d(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%dp%d_b1_bx(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%dp%d_b1_by(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n",
      i,j,i,j,  i,j,i,j, i,j,i,j
     );

    }
  }

  printf("/* BUBBLE ( 0 , 1 ) */\n\n");

  for (int i = 2; i <= p1 + 1; i++)
  {
    for (int j = 0; j <= p2; j++)
    {
     printf(
      "static double gradleg_quad_p%dp%d_b2_a(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%dp%d_b2_b(double x, double y)\n"
      "{\n"
      "  return l%d(x) * Legendre%d(y);\n"
      "}\n\n",
      i,j,  i,j,i,j
     );
     printf(
      "static double gradleg_quad_p%dp%d_b2_ax(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%dp%d_b2_ay(double x, double y)\n"
      "{\n"
      "  return 0.0;\n"
      "}\n\n"
      "static double gradleg_quad_p%dp%d_b2_bx(double x, double y)\n"
      "{\n"
      "  return dl%d(x) * Legendre%d(y);\n"
      "}\n\n"
      "static double gradleg_quad_p%dp%d_b2_by(double x, double y)\n"
      "{\n"
      "  return l%d(x) * Legendre%dx(y);\n"
      "}\n\n",
       i,j,i,j,  i,j,i,j, i,j,i,j
     );

    }
  }
}

int main(int argc, char* argv[])
{
  int i, j, k, l;
  printf("/* Whitney fns - constant tangential component */\n\n");
  whitney(0);

  printf("/* Edge fns - gradients of scalar lobatto edge functions */\n\n");
  for (i = 0; i <= 1; i++)
    for (j = 2; j <= 11; j++)
      edge_fn(i, j);
  for (j = 0; j <= 1; j++)
    for (i = 2; i <= 11; i++)
      edge_fn(i, j);


  bubble(10,10);

  printf("static Shapeset::shape_fn_t gradleg_quad_fn_a[] = \n{\n");      
  
  printf("  gradleg_quad_p0_e1_a_0, gradleg_quad_p0_e1_a_1, gradleg_quad_p0_e2_a, gradleg_quad_p0_e2_a, gradleg_quad_p0_e3_a_0, gradleg_quad_p0_e3_a_1, gradleg_quad_p0_e4_a, gradleg_quad_p0_e4_a, \n");

  for (j = 2; j <= 11; j++)
  {
    if ((j%2))
      printf("  gradleg_quad_l%d_l0_a_0, gradleg_quad_l%d_l0_a_1, gradleg_quad_l1_l%d_a_0, gradleg_quad_l1_l%d_a_1,  gradleg_quad_l%d_l1_a_0, gradleg_quad_l%d_l1_a_1, gradleg_quad_l0_l%d_a_0, gradleg_quad_l0_l%d_a_1, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradleg_quad_l%d_l0_a, gradleg_quad_l%d_l0_a, gradleg_quad_l1_l%d_a, gradleg_quad_l1_l%d_a, gradleg_quad_l%d_l1_a, gradleg_quad_l%d_l1_a, gradleg_quad_l0_l%d_a, gradleg_quad_l0_l%d_a, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }
  printf("\n");

  k = 88;
  int indices1[15][15];
  int indices2[15][15];
  for (int i = 0; i <= 11; i++)
    for (int j = 0; j <= 11; j++)
    {
      indices1[i][j] = 0;
      indices2[i][j] = 0;
    }
 
  for (int i = 0; i <= 10; i++)
    for (int j = 2; j <= 10 + 1; j++)
    {
      printf("  gradleg_quad_p%dp%d_b1_a, ", i,j);
      indices1[i][j-1] = k;
      k++;
    }
  for (int i = 2; i <= 10 + 1; i++)
    for (int j = 0; j <= 10; j++)
    {
      printf("  gradleg_quad_p%dp%d_b2_a, ", i,j);
      indices2[i-1][j] = k;
      k++;    
    } 
  printf("};\n\n");




  printf("static Shapeset::shape_fn_t gradleg_quad_fn_b[] = \n{\n");  

  printf("  gradleg_quad_p0_e1_b, gradleg_quad_p0_e1_b, gradleg_quad_p0_e2_b_0, gradleg_quad_p0_e2_b_1,  gradleg_quad_p0_e3_b, gradleg_quad_p0_e3_b, gradleg_quad_p0_e4_b_0, gradleg_quad_p0_e4_b_1, \n");
    
  for (j = 2; j <= 11; j++)
  {
    if (!(j%2))
      printf("  gradleg_quad_l%d_l0_b, gradleg_quad_l%d_l0_b, gradleg_quad_l1_l%d_b, gradleg_quad_l1_l%d_b, gradleg_quad_l%d_l1_b, gradleg_quad_l%d_l1_b, gradleg_quad_l0_l%d_b,  gradleg_quad_l0_l%d_b, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradleg_quad_l%d_l0_b_0, gradleg_quad_l%d_l0_b_1, gradleg_quad_l1_l%d_b_0, gradleg_quad_l1_l%d_b_1, gradleg_quad_l%d_l1_b_0, gradleg_quad_l%d_l1_b_1, gradleg_quad_l0_l%d_b_0, gradleg_quad_l0_l%d_b_1, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }  
  printf("\n");

  for (int i = 0; i <= 10; i++)
    for (int j = 2; j <= 10 + 1; j++)
      printf("  gradleg_quad_p%dp%d_b1_b, ", i,j);
  for (int i = 2; i <= 10 + 1; i++)
    for (int j = 0; j <= 10; j++)
      printf("  gradleg_quad_p%dp%d_b2_b, ", i,j);

  printf("};\n\n");


  printf("static Shapeset::shape_fn_t gradleg_quad_fn_ax[] = \n{\n");      
  
  printf("  gradleg_quad_p0_e1_ax_0, gradleg_quad_p0_e1_ax_1, gradleg_quad_p0_e2_ax, gradleg_quad_p0_e2_ax, gradleg_quad_p0_e3_ax_0, gradleg_quad_p0_e3_ax_1, gradleg_quad_p0_e4_ax, gradleg_quad_p0_e4_ax, \n");

  for (j = 2; j <= 11; j++)
  {
    if ((j%2))
      printf("  gradleg_quad_l%d_l0_ax_0, gradleg_quad_l%d_l0_ax_1, gradleg_quad_l1_l%d_ax_0, gradleg_quad_l1_l%d_ax_1,  gradleg_quad_l%d_l1_ax_0, gradleg_quad_l%d_l1_ax_1, gradleg_quad_l0_l%d_ax_0, gradleg_quad_l0_l%d_ax_1, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradleg_quad_l%d_l0_ax, gradleg_quad_l%d_l0_ax, gradleg_quad_l1_l%d_ax, gradleg_quad_l1_l%d_ax, gradleg_quad_l%d_l1_ax, gradleg_quad_l%d_l1_ax, gradleg_quad_l0_l%d_ax, gradleg_quad_l0_l%d_ax, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }
  printf("\n");
  for (int i = 0; i <= 10; i++)
    for (int j = 2; j <= 10 + 1; j++)
      printf("  gradleg_quad_p%dp%d_b1_ax, ", i,j);
  for (int i = 2; i <= 10 + 1; i++)
    for (int j = 0; j <= 10; j++)
      printf("  gradleg_quad_p%dp%d_b2_ax, ", i,j);

  printf("};\n\n");




  printf("static Shapeset::shape_fn_t gradleg_quad_fn_bx[] = \n{\n");  

  printf("  gradleg_quad_p0_e1_bx, gradleg_quad_p0_e1_bx, gradleg_quad_p0_e2_bx_0, gradleg_quad_p0_e2_bx_1,  gradleg_quad_p0_e3_bx, gradleg_quad_p0_e3_bx, gradleg_quad_p0_e4_bx_0, gradleg_quad_p0_e4_bx_1, \n");
    
  for (j = 2; j <= 11; j++)
  {
    if (!(j%2))
      printf("  gradleg_quad_l%d_l0_bx, gradleg_quad_l%d_l0_bx, gradleg_quad_l1_l%d_bx, gradleg_quad_l1_l%d_bx, gradleg_quad_l%d_l1_bx, gradleg_quad_l%d_l1_bx, gradleg_quad_l0_l%d_bx,  gradleg_quad_l0_l%d_bx, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradleg_quad_l%d_l0_bx_0, gradleg_quad_l%d_l0_bx_1, gradleg_quad_l1_l%d_bx_0, gradleg_quad_l1_l%d_bx_1, gradleg_quad_l%d_l1_bx_0, gradleg_quad_l%d_l1_bx_1, gradleg_quad_l0_l%d_bx_0, gradleg_quad_l0_l%d_bx_1, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }  
  printf("\n");
  for (int i = 0; i <= 10; i++)
    for (int j = 2; j <= 10 + 1; j++)
      printf("  gradleg_quad_p%dp%d_b1_bx, ", i,j);
  for (int i = 2; i <= 10 + 1; i++)
    for (int j = 0; j <= 10; j++)
      printf("  gradleg_quad_p%dp%d_b2_bx, ", i,j);

  printf("};\n\n");


  printf("static Shapeset::shape_fn_t gradleg_quad_fn_ay[] = \n{\n");      
  
  printf("  gradleg_quad_p0_e1_ay_0, gradleg_quad_p0_e1_ay_1, gradleg_quad_p0_e2_ay, gradleg_quad_p0_e2_ay, gradleg_quad_p0_e3_ay_0, gradleg_quad_p0_e3_ay_1, gradleg_quad_p0_e4_ay, gradleg_quad_p0_e4_ay, \n");

  for (j = 2; j <= 11; j++)
  {
    if ((j%2))
      printf("  gradleg_quad_l%d_l0_ay_0, gradleg_quad_l%d_l0_ay_1, gradleg_quad_l1_l%d_ay_0, gradleg_quad_l1_l%d_ay_1,  gradleg_quad_l%d_l1_ay_0, gradleg_quad_l%d_l1_ay_1, gradleg_quad_l0_l%d_ay_0, gradleg_quad_l0_l%d_ay_1, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradleg_quad_l%d_l0_ay, gradleg_quad_l%d_l0_ay, gradleg_quad_l1_l%d_ay, gradleg_quad_l1_l%d_ay, gradleg_quad_l%d_l1_ay, gradleg_quad_l%d_l1_ay, gradleg_quad_l0_l%d_ay, gradleg_quad_l0_l%d_ay, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }
  printf("\n");
   for (int i = 0; i <= 10; i++)
    for (int j = 2; j <= 10 + 1; j++)
      printf("  gradleg_quad_p%dp%d_b1_ay, ", i,j);
  for (int i = 2; i <= 10 + 1; i++)
    for (int j = 0; j <= 10; j++)
      printf("  gradleg_quad_p%dp%d_b2_ay, ", i,j);

  printf("};\n\n");




  printf("static Shapeset::shape_fn_t gradleg_quad_fn_by[] = \n{\n");  

  printf("  gradleg_quad_p0_e1_by, gradleg_quad_p0_e1_by, gradleg_quad_p0_e2_by_0, gradleg_quad_p0_e2_by_1,  gradleg_quad_p0_e3_by, gradleg_quad_p0_e3_by, gradleg_quad_p0_e4_by_0, gradleg_quad_p0_e4_by_1, \n");
    
  for (j = 2; j <= 11; j++)
  {
    if (!(j%2))
      printf("  gradleg_quad_l%d_l0_by, gradleg_quad_l%d_l0_by, gradleg_quad_l1_l%d_by, gradleg_quad_l1_l%d_by, gradleg_quad_l%d_l1_by, gradleg_quad_l%d_l1_by, gradleg_quad_l0_l%d_by,  gradleg_quad_l0_l%d_by, ", j,j,j,j,j,j,j,j);
    else
      printf("  gradleg_quad_l%d_l0_by_0, gradleg_quad_l%d_l0_by_1, gradleg_quad_l1_l%d_by_0, gradleg_quad_l1_l%d_by_1, gradleg_quad_l%d_l1_by_0, gradleg_quad_l%d_l1_by_1, gradleg_quad_l0_l%d_by_0, gradleg_quad_l0_l%d_by_1, ", j,j,j,j,j,j,j,j);
    printf("\n");
  }  
  printf("\n");
  int s = 88;
  for (int i = 0; i <= 10; i++)
    for (int j = 2; j <= 10 + 1; j++)
    { 
      printf("  gradleg_quad_p%dp%d_b1_by, ", i,j);
      s++;
    }
  for (int i = 2; i <= 10 + 1; i++)
    for (int j = 0; j <= 10; j++)
      printf("  gradleg_quad_p%dp%d_b2_by, ", i,j);


  printf("};\n\n");

////////////////////////////////////////////////////////////////////////////

 
  for (int i = 0; i <= 10; i++)
    for (int j = 0; j <= 10; j++)
    {
      if ((i != 0) || (j!= 0))
      {
        printf("static int qb_%d_%d[] = { ", i, j);
        for (int k = 0; k <= i; k++)
          for (int l = 0; l <= j; l++)
          {
            if (indices1[k][l] != 0) printf("%d,", indices1[k][l]);
            if (indices2[k][l] != 0) printf("%d,", indices2[k][l]);
          } 
        printf("};\n  ");      
      }
    }



  printf("\n\n");

  printf("#define NULL16 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,\n\n");
  printf("static int* gradleg_quad_bubble_indices[] =\n{\n");
  //printf("  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL, NULL16\n");
  for (i = 0; i < 1; i++)
  {
    printf("  NULL, ");
    for (j = 1; j <= 10; j++)
      printf("qb_%d_%d, ", i, j);
    printf(" NULL, NULL, NULL, NULL, NULL, NULL16\n");
  }
  for (i = 1; i <= 10; i++)
  {
    printf("  ");
    for (j = 0; j <= 10; j++)
      printf("qb_%d_%d, ", i, j);
    printf(" NULL, NULL, NULL, NULL, NULL, NULL16\n");
  }
  printf("};\n\n");
 
  printf("#define zero16  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,\n\n");
  printf("static int gradleg_quad_bubble_count[] =\n{\n");
  //printf("  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, zero16\n");
  for (i = 0; i < 1; i++)
  {  
    printf("  0,  ");
    for (j = 1; j <= 10; j++)
    {
      printf("%d,  ", j);
    }
    printf("0,  0,  0,  0,  0, zero16 \n");
  }
  for (i = 1; i <= 10; i++)
  {  
    printf("  ");
    for (j = 0; j <= 10; j++)
    {
      printf("%d, ", (j+1)*i +(i+1)*j);
    }
    printf("0,  0,  0,  0,  0, zero16 \n");
  }
  printf("};\n\n");

  printf("static int gradleg_quad_vertex_indices[4] ={-1, -1, -1, -1};\n\n");

  for (int edge = 0; edge < 4; edge++)
  {
    printf("static int gradleg_quad_edge_indices_%d[] = { ", edge);
    k = 0;
    for (i = 0; i <= 10; i++)
    {
      printf("%d, %d, ", 2*edge + k, 2*edge + k +1);
      k += 8;
    }
    printf("};\n\n");
  }
  printf("static int* gradleg_quad_edge_indices[4] =\n{\n");
  for (int i = 0; i < 4; i++)
    printf("  gradleg_quad_edge_indices_%d,\n",i);
  printf("};\n\n");

  printf("#define oo make_quad_order\n\n");

  printf("static int gradleg_quad_index_to_order[] = \n{\n");      
  
  printf("  oo(0,1), oo(0,1), oo(1,0), oo(1,0), oo(0,1), oo(0,1), oo(1,0), oo(1,0),\n");

  for (j = 2; j <= 11; j++)
  {
    printf("  oo(%d,1), oo(%d,1), oo(1,%d), oo(1,%d), oo(%d,1), oo(%d,1), oo(1,%d), oo(1,%d),\n", j,j,j,j,j,j,j,j);
  }
  printf("\n");
  for (int i = 0; i <= 10; i++)
    for (int j = 2; j <= 10 + 1; j++)
      printf("  oo(%d,%d),",i,j);
  for (int i = 2; i <= 10 + 1; i++)
    for (int j = 0; j <= 10; j++)
      printf("  oo(%d,%d),",i,j);

  printf("\n};\n\n");


}

