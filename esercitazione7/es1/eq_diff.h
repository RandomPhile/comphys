#ifndef eq_diff_h
#define eq_diff_h
#include "funzioni_utili.h"

using namespace std;

extern int N_mol;

void primi_vicini(double *r, double *r_prim_vic, double L, int N_mol, int j) {
    for (int i = 0; i < 3 * N_mol; i += 3) {
        if (i / 3 == j) { //se la particella è se stessa impongo che la distanza dal suo primo vicino sia 0
            r_prim_vic[i] = 0;
            r_prim_vic[i + 1] = 0;
            r_prim_vic[i + 2] = 0;
        }
        else { //trova il primo vicino
            r_prim_vic[i]     = r[i]     - r[3 * j]     + L * rint((r[3 * j]     - r[i])     / L);
            r_prim_vic[i + 1] = r[i + 1] - r[3 * j + 1] + L * rint((r[3 * j + 1] - r[i + 1]) / L);
            r_prim_vic[i + 2] = r[i + 2] - r[3 * j + 2] + L * rint((r[3 * j + 2] - r[i + 2]) / L);
        }
    }
}

double VLJ(double r, double L) {//trova il valore del potenziale che ho da una particella
    if (r < L / 2 && r != 0) {
        double Vrc = 4 * (pow1(2 / L, 12) - pow1(2 / L, 6));
        return 4 * (pow1(1 / r, 12) - pow1(1 / r, 6)) - Vrc;
    } else {
        return 0;
    }
}

void fLJ(double *r, double L, long double *F, int i) { //arg[0]=eps, arg[1]=sigma, arg[2]= dimensione interazione, arg[3]=dimensione scatola
    setta_matr(F, 0, 1);//dim F = 3
    
    double r_pv[3 * N_mol];
    setta_matr(r_pv, 0, N_mol);
    
    primi_vicini(r, r_pv, L, N_mol, i);//modifica r_pv

    double mod_r_pv[N_mol];

    for (int j = 0; j < 3 * N_mol; j += 3) {
        mod_r_pv[j / 3] = sqrt(r_pv[j] * r_pv[j] + r_pv[j + 1] * r_pv[j + 1] + r_pv[j + 2] * r_pv[j + 2]);
    }


    for (int j = 0; j < 3 * N_mol; j += 3) {
        if ((r_pv[j] < L / 2) && (r_pv[j + 1] < L / 2) && (r_pv[j + 2] < L / 2))  { //se è abbastanza vicina <L/2
            if (j / 3 != i) {
                long double r7 = pow1((long double) mod_r_pv[j / 3], 7);
                long double cost = 24 * (2 / (r7 * r7) - 1 / (r7 * mod_r_pv[j / 3]));
                F[0] += (long double) r_pv[j] * cost;
                F[1] += (long double) r_pv[j + 1] * cost;
                F[2] += (long double) r_pv[j + 2] * cost;
            }
        }
    }
}

int vel_verlet(double t, double *r, double *v, double dt, int j,
               double *a_prev,
               void (*fLJ)(double*, double , long double*, int), double L) {//per potenziale generico
    long double F[3];
    double temp[3];
    for (int i = 3 * j; i < 3 * j + 3; ++i) {
        temp[i - 3 * j] = a_prev[i];
        r[i] = r[i] + v[i] * dt + a_prev[i] * dt * dt / 2;
    }

    fLJ(r, L, F, j);//aggiorna F[3]

    for (int i = 3 * j; i < 3 * j + 3; ++i) {
        a_prev[i] = F[i - 3 * j];
        v[i] += dt * (a_prev[i] + temp[i - 3 * j]) / 2;
    }
    return 0;
}

#endif /* eq_diff_h */