#include <stdio.h>
#include <math.h>

int eulero_exp(double t, double h, double *x, double *y, double(*f)(double, double, double, double*),
               double *fargs, double(*g)(double, double, double, double* ), double *gargs) {
    
    double x_n = *x;           
   *x = *x + h * f(t, *x, *y, fargs);
   *y = *y + h * g(t, x_n, *y, gargs);               
               
               


   return 0;
}

int runge_kutta (double t, double h, double *x, double *y, double(*f)(double, double, double, double*),
               double *fargs, double(*g)(double, double, double, double* ), double *gargs) {
               
               
   double k1, k2, k3, k4;
   double l1, l2, l3, l4;
   
   k1 = h * f(t, *x, *y, fargs);
   l1 = h * g(t, *x, *y, gargs);
   
   k2 = h * f(t + h/2, *x + k1/2, *y + l1/2, fargs);
   l2 = h * g(t + h/2, *x + k1/2, *y + l1/2, gargs);
   
   k3 = h * f(t + h/2, *x + k2/2, *y + l2/2, fargs);
   l3 = h * g(t + h/2, *x + k2/2, *y + l2/2, gargs);
   
   k4 = h * f(t + h, *x + k3, *y + l3, fargs);
   l4 = h * g(t + h, *x + k3, *y + l3, gargs);
   
   *x = *x + ( k1 + 2*k2 + 2*k3 + k4 )/6;
   *y = *y + ( l1 + 2*l2 + 2*l3 + l4 )/6;            
               
   return 0;               
               
}            
               
               
               
               
               
               
               
               
               
               
               
               








