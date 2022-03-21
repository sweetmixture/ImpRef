#include<math.h>

// Rosenbrock

double RosenBrock( const double a, const double b, const double x, const double y )
{
    return pow((a-x),2.) + b*pow((y-x*x),2.);
}

void D_RosenBrock( const double a, const double b, const double x, const double y, double* ref )
{
    double dx = -2.*(a-x) - 4.*x*b*(y-x*x);
    double dy = 2.*b*(y-x*x);
    
    ref[0] = dx;
    ref[1] = dy;
    
    return;
}
