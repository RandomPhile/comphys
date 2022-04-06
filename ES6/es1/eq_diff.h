#ifndef eq_diff_h
#define eq_diff_h

using namespace std;

extern double lam;

int eulero_exp(double t, double *x, double *y, double h,
               double (*f)(double, double, double, double*),
               double *fargs,
               double (*g)(double, double, double, double*),
               double *gargs) {
    double temp = *x;
    *x += h * f(t, *x, *y, NULL);
    *y += h * g(t, temp, *y, NULL);
    return 0;
}
int eulero_cromer(double t, double *x, double *y, double h,
                  double (*f)(double, double, double, double*),
                  double *fargs,
                  double (*g)(double, double, double, double*),
                  double *gargs) {
    *y += h * g(t, *x, *y, NULL);
    *x += h * f(t, *x, *y, NULL);
    return 0;
}
int eulero_imp(double t, double *x, double *y, double h,
               double (*f)(double, double, double, double*),
               double *fargs,
               double (*g)(double, double, double, double*),
               double *gargs) {
    double temp = *x;
    *x = ((*x) + h * f(t, *x, *y, NULL)) / (1 + pow(lam * h, 2));
    *y = ((*y) + pow(lam, 2) * h * g(t, temp, *y, NULL)) / (1 + pow(lam * h, 2));
    return 0;
}
int runge_kutta(double t, double *x, double *y, double h,
                double (*f)(double, double, double, double*),
                double *fargs,
                double (*g)(double, double, double, double*),
                double *gargs) {
    double k1 = h * f(t, *x, *y, NULL);
    double l1 = h * g(t, *x, *y, NULL);

    double k2 = h * f(t + h / 2, *x + k1 / 2, *y + l1 / 2, NULL);
    double l2 = h * g(t + h / 2, *x + k1 / 2, *y + l1 / 2, NULL);

    double k3 = h * f(t + h / 2, *x + k2 / 2, *y + l2 / 2, NULL);
    double l3 = h * g(t + h / 2, *x + k2 / 2, *y + l2 / 2, NULL);

    double k4 = h * f(t + h, *x + k3, *y + l3, NULL);
    double l4 = h * g(t + h, *x + k3, *y + l3, NULL);

    *x += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    *y += (l1 + 2 * l2 + 2 * l3 + l4) / 6;
    return 0;
}

int vel_verlet(double t, double *x, double *v, double dt,
               double *a_prev,
               double (*a)(double, double*),
               double *args) {
    double temp = *a_prev;
    *x = *x + *v * dt + a(*x, args) * pow(dt, 2) / 2;
    *a_prev = - *x;
    *v = *v + dt * (*a_prev + temp) / 2;
    return 0;
}
#endif /* eq_diff_h */