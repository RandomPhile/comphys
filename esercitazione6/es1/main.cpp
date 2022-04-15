#include <iostream>
#include <fstream>
#include <cmath>
#include "eq_diff.h"

using namespace std;
ofstream dati;

//putatore NULL
double f(double t, double x, double y, double *fargs);
double g(double t, double x, double y, double *gargs);
double a(double *r, int indice, double *args);

double lam;

int main() {
    dati.open("dati.dat");
    double t , h, x = 1, y = 0, E,K, T;
    double t0 = 0, t1 = 25;
    int N_mol = 10000;
    int N = 1000;
    h = (t1 - t0) / ( (double) N);

    double r[3 * N_mol], v[3 * N_mol];
    //resetta_matr(r, 0, N_mol);
    gauss_distr(r, 1, N_mol);
    gauss_distr(v, 1, N_mol);
    set_vcm0(v,N_mol);

    double a_prev[3 * N_mol];
    resetta_matr(a_prev, r, N_mol);
    //inizializzo accelerazioni come -r

    double r_mod, v_mod;

    E = 0; K = 0;
    for (int n = 0; n < N; ++n) {//passi temporali
        E = E*(n+1);//energia totale @t
        K = K*(n+1);//energia cinetica @t
        T = 0;//temp @t
        for (int j = 0; j < N_mol; ++j) {//particelle
            t = t0 + n * h;
            if (vel_verlet(t, r, v, h, j, a_prev, a, NULL)) {printf("ERRORE");}

            r_mod = sqrt(r[j] * r[j] + r[j + 1] * r[j + 1] + r[j + 2] * r[j + 2]);
            v_mod = sqrt(v[j] * v[j] + v[j + 1] * v[j + 1] + v[j + 2] * v[j + 2]);
            K += 0.5 * v_mod * v_mod;
            E += 0.5 * r_mod * r_mod + 0.5 * v_mod * v_mod;
            //T += v_mod * v_mod * 0.5 / N_mol;
        }
        E = E/(n+2);
        K = K/(n+2);
        T = 2.0*K/(3.0*N_mol);
        dati << t  << "\t" << E << "\t" << K << "\t" << T << endl;
    }
    dati.close();
    return 0;
}

double f(double t, double x, double y, double *fargs) {
    return y;
}
double g(double t, double x, double y, double *gargs) {
    return -x;
}
double a(double *r, int indice, double *args) {
    return -r[indice];
}

