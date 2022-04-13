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

//int vel_verlet(double t, double r[][3], double v[][3], double dt, int N_Mol,
//               double f_prev[][3],
//               double (*f)(double*, int, double*),
//               double *fargs) {
//    
//    for (int j = 0; j < 3; ++j) {
//        double temp[N_Mol];
//        for (int i=0;i<N_Mol;i++) {
//            temp[i] = f_prev[i][j];
//            r[i][j] = r[i][j] + v[i][j] * dt + f(r,i,fargs)* pow(dt, 2) / 2;
//            f_prev[i][j] = -r[i][j];
//            v[i][j] += dt*(f_prev[i][j] + temp[i]) / 2;
//        }
//    }
//    
//    return 0;
//}

#endif /* eq_diff_h */
