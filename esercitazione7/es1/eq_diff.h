#ifndef eq_diff_h
#define eq_diff_h
using namespace std;

int vel_verlet(double t, double *r, double *v, double dt, int j, int dim,
               double *a_prev,
               double (*a)(double*, int, double*),
               double *args) {
    if(dim==3){
        for (int i = 3 * j; i < 3 * j + 3; ++i) {
            double temp = a_prev[i];
            r[i] = r[i] + v[i] * dt + a(r, i, args) * pow(dt, 2) / 2;
            a_prev[i] = -r[i];
            v[i] += dt * (a_prev[i] + temp) / 2;
        }
    }
    else{
        for (int i = 3 * j; i < 3 * j + 3; ++i) {
            double temp = a_prev[i];
            r[i] = r[i] + v[i] * dt + a(r, i, args) * pow(dt, 2) / 2;
            a_prev[i] = -r[i];
            v[i] += dt * (a_prev[i] + temp) / 2;
        }
    }
    return 0;
}
double pow1(double base, int esp){
    double ris=1;
    for(int i=0; i<esp; i++){
        ris*=base;
    }
    return ris;
}
int vel_verlet(double t, double *r, double *v, double dt, int j, int dim,
               double *a_prev,
               double (*f)(double*, double*, int),
               double *args) {//per potenziale generico
    if(dim==3){
        for (int i = 3 * j; i < 3 * j + 3; ++i) {
            double temp = a_prev[i];
            r[i] = r[i] + v[i] * dt + f(r, args, i) * pow(dt, 2) / 2;
            a_prev[i] = f(r, args, i);
            v[i] += dt * (a_prev[i] + temp) / 2;
        }
    }
    else{
        for (int i = j; i < j + 1; ++i) {
            double temp = a_prev[i];
            r[i] = r[i] + v[i] * dt + f(r, args, i) * pow(dt, 2) / 2;
            a_prev[i] = f(r, args, i);
            v[i] += dt * (a_prev[i] + temp) / 2;
        }
    }
    return 0;
}
#endif /* eq_diff_h */
