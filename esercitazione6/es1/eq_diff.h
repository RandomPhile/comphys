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

int vel_verlet(double t, double *r, double *v, double dt, int j,
               double *a_prev,
               double (*a)(double*, int, double*),
               double *args) {
    for (int i = 3 * j; i < 3 * j + 3; ++i) {
        double temp = a_prev[i];
        r[i] = r[i] + v[i] * dt + a(r, i, args) * pow(dt, 2) / 2;
        a_prev[i] = -r[i];
        v[i] += dt * (a_prev[i] + temp) / 2;
    }
    return 0;
}

void gauss_distr(double *v, double sigma, int N_Mol) {
    for (int i = 0; i < 3 * N_Mol; i += 2) {
        double x1, x2;
        x1 = rand() / ((double)RAND_MAX + 1.0);
        x2 = rand() / ((double)RAND_MAX + 1.0);
        v[i] = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
        v[i + 1] = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1);
    }
}

void resetta_matr(double *r, int val, int N_Mol) {
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = val;
    }
}
void resetta_matr(double *r, double *val, int N_Mol) {
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = -val[j];
    }
}
void set_vcm0(double *v, int N_Mol){
    double v_cm[]={0,0,0};
    for(int j=0;j<3;j++){
        for(int i=0;i<N_Mol;i++){
            v_cm[j]+=v[i+3*j];
        }
    }
    for(int j=0;j<3;j++){
        v_cm[j]/=N_Mol;
    }
    for(int j=0;j<3;j++){
        for(int i=0;i<N_Mol;i++){
            v[i+3*j]-=v_cm[j];
        }
    }
}
#endif /* eq_diff_h */
