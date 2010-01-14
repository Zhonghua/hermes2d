#include "numerical_flux.h"

double matrix_R(int i, int j, double w0, double w1, double w3, double w4)
{
    double rho = w0;
    double u = w1/w0;
    double w = w3/w0;
    double E = w4;
    double v2 = u*u+w*w;
    double kappa=1.4;
    double p = (kappa-1)*(E - rho*v2/2);
    double c=sqrt(kappa*p/rho);
    if (i == 0 && j == 0)
        return 1;
    else if (i == 0 && j == 1)
        return 1;
    else if (i == 0 && j == 2)
        return 1;
    else if (i == 0 && j == 3)
        return 1;

    else if (i == 1 && j == 0)
        return u-c;
    else if (i == 1 && j == 1)
        return u;
    else if (i == 1 && j == 2)
        return u;
    else if (i == 1 && j == 3)
        return u+c;

    else if (i == 2 && j == 0)
        return w;
    else if (i == 2 && j == 1)
        return w;
    else if (i == 2 && j == 2)
        return w-c;
    else if (i == 2 && j == 3)
        return w;

    else if (i == 3 && j == 0)
        return v2/2 + c*c/(kappa-1) - u*c;
    else if (i == 3 && j == 1)
        return v2/2;
    else if (i == 3 && j == 2)
        return v2/2 - w*c;
    else if (i == 3 && j == 3)
        return v2/2 + c*c/(kappa-1) + u*c;

    printf("i=%d, j=%d;\n", i, j);
    error("Invalid index.");
}

double matrix_R_inv(int i, int j, double w0, double w1, double w3, double w4)
{
    double rho = w0;
    double u = w1/w0;
    double w = w3/w0;
    double E = w4;
    double v2 = u*u+w*w;
    double kappa=1.4;
    double p = (kappa-1)*(E - rho*v2/2);
    double c=sqrt(kappa*p/rho);
    if (i == 0 && j == 0)
        return ((kappa-1)*v2/2 + u*c)/2;
    else if (i == 0 && j == 1)
        return -(c+u*(kappa-1))/2;
    else if (i == 0 && j == 2)
        return -w*(kappa-1)/2;
    else if (i == 0 && j == 3)
        return (kappa-1)/2;

    else if (i == 1 && j == 0)
        return c*c-c*w-(kappa-1)*v2/2;
    else if (i == 1 && j == 1)
        return u*(kappa-1);
    else if (i == 1 && j == 2)
        return c+w*(kappa-1);
    else if (i == 1 && j == 3)
        return 1-kappa;

    else if (i == 2 && j == 0)
        return w*c;
    else if (i == 2 && j == 1)
        return 0;
    else if (i == 2 && j == 2)
        return -c;
    else if (i == 2 && j == 3)
        return 0;

    else if (i == 3 && j == 0)
        return ((kappa-1)*v2/2 - u*c)/2;
    else if (i == 3 && j == 1)
        return (c-u*(kappa-1))/2;
    else if (i == 3 && j == 2)
        return -w*(kappa-1)/2;
    else if (i == 3 && j == 3)
        return (kappa-1)/2;

    printf("i=%d, j=%d;\n", i, j);
    error("Invalid index.");
}
