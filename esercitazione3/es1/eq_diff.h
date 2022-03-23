//
//  eq_diff.h
//  es1
//
//  Created by Mattia Lupi on 16/03/22.
//

#ifndef eq_diff_h
#define eq_diff_h

using namespace std;

extern double lam;

int eulero_exp(double t, double *x, double *y, double h,
               double (*f)(double,double,double,double*),
               double *fargs,
               double (*g)(double,double,double,double*),
               double *gargs){
    double temp= *x;
    *x += h*f(t,*x,*y,NULL);
    *y += h*g(t,temp,*y,NULL);
    return 0;
}
int eulero_cromer(double t, double *x, double *y, double h,
               double (*f)(double,double,double,double*),
               double *fargs,
               double (*g)(double,double,double,double*),
               double *gargs){
    *y+=h*g(t,*x,*y,NULL);
    *x+=h*f(t,*x,*y,NULL);
    return 0;
}
int eulero_imp(double t, double *x, double *y, double h,
               double (*f)(double,double,double,double*),
               double *fargs,
               double (*g)(double,double,double,double*),
               double *gargs){
    double temp= *x;
    *x=((*x)+h*f(t,*x,*y,NULL))/(1+pow(lam*h,2));
    *y=((*y)+pow(lam,2)*h*g(t,temp,*y,NULL))/(1+pow(lam*h,2));
    return 0;
}

#endif /* eq_diff_h */
