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

void gauss_distr(double *v, double sigma, int N_Mol) {
    for (int i = 0; i < 3 * N_Mol; i += 2) {
        double x1, x2;
        x1 = rand() / ((double)RAND_MAX + 1.0);
        x2 = rand() / ((double)RAND_MAX + 1.0);
        v[i] = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
        v[i + 1] = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1);
    }
}

void setta_matr(double *r, int val, int N_Mol) {
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = val;
    }
}
void compila_matr(double *r,int N_Mol,
                  double (*f)(double*, double*, double, double),
                  double *fargs, double arg) {
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = f(r,fargs,arg,j);
    }
}
void copia_vett(double *r, int N_Mol,double *v) {
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = v[j];
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
