#include <iostream>
#include <fstream>
#include <cmath>
#include "eq_diff.h"
#include "funzioni_utili.h"

using namespace std;
ofstream dati;


double f(double t, double x, double y, double *fargs);
double g(double t, double x, double y, double *gargs);
double a(double *r, int indice, double *args);

int N_mol = pow1(2, 3);

int main() {
    dati.open("dati.dat");
    double t = 0, h, E, K, T,V;
    double t0 = 0, t1 = 0.01;

    int N = 10;
    h = (t1 - t0) / ( (double) N);
    double rho = 1e-2; //fisso la densita del campione da studiare
    double eps, sigma, L, distanza_interaz;
    L = cbrt(N_mol) / rho;
    distanza_interaz = L / 2;

    double var_ad[] = {eps, sigma, distanza_interaz, L};

    double r[3 * N_mol], v[3 * N_mol], r_prim_vic[3 * N_mol];
    double r_mod, v_mod;

    r_mod = sqrt(r[0] * r[0] + r[0 + 1] * r[0 + 1] + r[0 + 2] * r[0 + 2]);
    v_mod = sqrt(v[0] * v[0] + v[0 + 1] * v[0 + 1] + v[0 + 2] * v[0 + 2]);

    crea_reticolo(N_mol, L, r);
    
    gauss_distr(v, 1, N_mol);
    set_vcm0(v, N_mol);

    double a_prev[3 * N_mol];
    compila_matr(r, a_prev, N_mol, fLJ, var_ad);//inizializzo accelerazioni come date dal potenziale
    
    E = 0; K = 0; V = 0;
    for (int n = 0; n < N; ++n) {//passi temporali
        //E = E*(n+1);//energia totale @t
        K = K * (n + 1); //energia cinetica @t
        V = V * (n + 1); //energia potenziale @t
        T = 0;//temp @t

        r_mod = sqrt(r[0] * r[0] + r[0 + 1] * r[0 + 1] + r[0 + 2] * r[0 + 2]);
        v_mod = sqrt(v[0] * v[0] + v[0 + 1] * v[0 + 1] + v[0 + 2] * v[0 + 2]);

        for (int j = 0; j < N_mol; ++j) {//particelle
            t = t0 + n * h;
            cout<<"tempo= "<<t<<endl;
            if (vel_verlet(t, r, v, h, j, 3, a_prev, fLJ, var_ad)) {printf("ERRORE");}
            
            r_mod = sqrt(r[j] * r[j] + r[j + 1] * r[j + 1] + r[j + 2] * r[j + 2]);
            v_mod = sqrt(v[j] * v[j] + v[j + 1] * v[j + 1] + v[j + 2] * v[j + 2]);

            cond_bordo(r, j, L);
            for (int i = 0; i < j; ++i) {//attenzione che è usata la distanza effettiva invece della distanza della prima copia vicina 
                double x_ij=r[3 * i] - r[3 * j];
                double y_ij=r[3 * i + 1] - r[3 * j + 1];
                double z_ij=r[3 * i + 2] - r[3 * j + 2];
                V += VLJ(sqrt(x_ij * x_ij + y_ij * y_ij + z_ij * z_ij), var_ad);
            }
            K += 0.5 * v_mod * v_mod;

        }
        // E = E/(n+2);
        K = K / (n + 2);
        V = V / (n + 2);

        E = K + V;
        T = 2.0 * K / (3.0 * N_mol);
        //dati << t  << "\t" << E << "\t" << K << "\t" << V << "\t" << T << endl;
        //cout << t  << "\t" << E << "\t" << K << "\t" << V << "\t" << T << endl;
    }
    dati<<"\n"<<endl;
    stampa_reticolo(N_mol, L, r);
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

