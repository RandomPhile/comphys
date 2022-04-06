#ifndef eq_diff_h
#define eq_diff_h
using namespace std;

int vel_verlet(double t, double *r, double *v, double dt,
               double *a_prev,
               double (*a)(double*, int, double*),
               double *args) {
    
    for (int i = 0; i < 3; ++i) {
        double temp = a_prev[i];
        r[i] = r[i] + v[i] * dt + a(r,i,args)* pow(dt, 2) / 2;
        a_prev[i] = -r[i];
        v[i] += dt*(a_prev[i] + temp) / 2;
    }
    return 0;
}
#endif /* eq_diff_h */