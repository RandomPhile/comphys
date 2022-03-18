//
//  eq_diff.h
//  es1
//
//  Created by Mattia Lupi on 16/03/22.
//

#ifndef eq_diff_h
#define eq_diff_h

using namespace std;


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
    
    *y+=h*g(t,*x,*y,NULL);
    *x+=h*f(t,*x,*y,NULL);
    return 0;
}

#endif /* eq_diff_h */
