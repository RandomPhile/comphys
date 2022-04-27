#ifndef eq_diff_h
#define eq_diff_h
#include "funzioni_utili.h"

using namespace std;

extern int N_mol;

void primi_vicini(double *r, double *r_prim_vic, double dim_scatola, int N_mol, int j) {
    
    for (int i = 0; i < 3 * N_mol; i += 3) {
        if(i/3==j){//se la particella è se stessa impongo che la distanza dal suo primo vicino sia 0
            r_prim_vic[i] = 0;
            r_prim_vic[i + 1]=0;
            r_prim_vic[i + 2]=0;
        }
        else{//trova il primo vicino
            r_prim_vic[i]     = r[i]     - r[3*j]     + dim_scatola * rint((r[3*j]     - r[i])     / dim_scatola);
            r_prim_vic[i + 1] = r[i + 1] - r[3*j + 1] + dim_scatola * rint((r[3*j + 1] - r[i + 1]) / dim_scatola);
            r_prim_vic[i + 2] = r[i + 2] - r[3*j + 2] + dim_scatola * rint((r[3*j + 2] - r[i + 2]) / dim_scatola);
        }
    }
}

double VLJ(double r, double *args) {//trova il valore del potenziale che ho da una particella
    if (r<args[2] && r!=0) {
        double Vrc = 4 * args[0] * (pow1(args[1] / args[2], 12) - pow1(args[1] / args[2], 6));
        return 4 * args[0] * (pow1(args[1] / r, 12) - pow1(args[1] / r, 6)) - Vrc;
    } else {
        return 0;
    }
}

void fLJ(double *r, double *args, long double *F, int i) { //arg[0]=eps, arg[1]=sigma, arg[2]= r_c dimensione interazione, arg[3]=dimensione scatola
    setta_matr(F, 0, 1);
    
    double r_pv[3 * N_mol];
    setta_matr(r_pv, 0, N_mol);
    primi_vicini(r, r_pv, args[3], N_mol, i);//modifica r_pv

    double mod_r_pv[N_mol];
    
    for (int j = 0; j < 3 * N_mol; j += 3) {
        mod_r_pv[j / 3] = sqrt(r_pv[j] * r_pv[j] + r_pv[j + 1] * r_pv[j + 1] + r_pv[j + 2] * r_pv[j + 2]);
    }
    

    for (int j = 0; j < 3 * N_mol; j += 3) {
        if ((r_pv[j] < args[2]) && (r_pv[j + 1] < args[2]) && (r_pv[j + 2] < args[2]))  { //se è abbastanza vicina <L/2
            if (j/3 != i) {
                long double r7 = pow1((long double) mod_r_pv[j/3], 7);
                long double sigma6 = pow1(args[1], 6);
                long double cost = 24 * args[0] * (2 * sigma6 * sigma6 / (r7 * r7) - sigma6 / (r7 * mod_r_pv[j/3]));
                F[0] += (long double) r_pv[j] * cost;
                F[1] += (long double) r_pv[j + 1] * cost;
                F[2] += (long double) r_pv[j + 2] * cost;
            }
        }
    }
}
int vel_verlet(double t, double *r, double *v, double dt, int j, int dim,
               double *a_prev,
               double (*a)(double*, int, double*),
               double *args) {//per oscillatori
    if (dim == 3) {
        for (int i = 3 * j; i < 3 * j + 3; ++i) {
            double temp = a_prev[i];
            r[i] = r[i] + v[i] * dt + a(r, i, args) * pow(dt, 2) / 2;
            a_prev[i] = -r[i];
            v[i] += dt * (a_prev[i] + temp) / 2;
        }
    }
    else {
        int i = j;
        double temp = a_prev[i];
        r[i] = r[i] + v[i] * dt + a(r, i, args) * pow(dt, 2) / 2;
        a_prev[i] = -r[i];
        v[i] += dt * (a_prev[i] + temp) / 2;
    }
    return 0;
}

int vel_verlet(double t, double *r, double *v, double dt, int j,
               double *a_prev,
               void (*fLJ)(double*, double*, long double*, int),
               double *args) {//per potenziale generico
    long double F[3];
    double temp[3];
    for (int i = 3 * j; i < 3 * j + 3; ++i) {
        temp[i - 3 * j] = a_prev[i];
        r[i] = r[i] + v[i] * dt + a_prev[i] * dt * dt / 2;
    }
    
    fLJ(r, args, F, j);//aggiorna F[3]
    
    for (int i = 3 * j; i < 3 * j + 3; ++i) {
        a_prev[i] = F[i - 3 * j];
        v[i] += dt * (a_prev[i] + temp[i - 3 * j]) / 2;
    }
    return 0;
}

#endif /* eq_diff_h */