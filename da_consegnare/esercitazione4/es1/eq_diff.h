#ifndef eq_diff_h
#define eq_diff_h
using namespace std;

int runge_kutta(double t, double *x, double *y, double h,
                double (*f)(double, double, double, double*),
                double *fargs,
                double (*g)(double, double, double, double*),
                double *gargs) {
    double k1 = h * f(t, *x, *y, fargs);
    double l1 = h * g(t, *x, *y, gargs);
    
    double k2 = h * f(t + h/2, *x + k1/2, *y + l1/2, fargs);
    double l2 = h * g(t + h/2, *x + k1/2, *y + l1/2, gargs);
    
    double k3 = h * f(t + h/2, *x + k2/2, *y + l2/2, fargs);
    double l3 = h * g(t + h/2, *x + k2/2, *y + l2/2, gargs);
    
    double k4 = h * f(t + h, *x + k3, *y + l3, fargs);
    double l4 = h * g(t + h, *x + k3, *y + l3, gargs);
    
    *x += (k1 + 2*k2 + 2*k3 + k4)/6;
    *y += (l1 + 2*l2 + 2*l3 + l4)/6;
    return 0;
}
int runge_kutta_r(double t, double *x, double *y, double h, double err,
                double (*f)(double, double, double, double*),
                double *fargs,
                double (*g)(double, double, double, double*),
                double *gargs) {
    double delta=pow(10,8);
    while(delta>err){
        double x1=*x, x2=*x, y1=*y, y2=*y;
        runge_kutta(t, x1, y1, h, f, fargs, g, gargs);
        runge_kutta(t, x2, y2, h/2, f, fargs, g, gargs);
        delta=fabs(x2-x1);
        if(delta!=0){
            h=min(0.9*(h*pow(err/delta,0.2)*4*h));
        }
        else{
            h*=4;
        }
        if(h<pow(10,-8)){
            break;
        }
    }
    
    
    return 0;
}

#endif /* eq_diff_h */
