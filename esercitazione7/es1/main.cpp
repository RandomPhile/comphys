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
double fLJ(double *r, double *args, int indice);


int main() {
    dati.open("dati.dat");
    double t=0, h, E, K, T;
    double t0 = 0, t1 = 25;
    int N_mol = pow1(10,3);
    int N = 1000;
    h = (t1 - t0) / ( (double) N);
    double rho=1e-2; //fisso la densita del campione da studiare
    double eps, sigma, L, distanza_interaz;
    L=rho*N_mol;
    
    double var_ad[]={eps, sigma, distanza_interaz};

    double r[3 * N_mol], v[3 * N_mol];
    
    crea_reticolo(N_mol, L, r);

    gauss_distr(v, 1, N_mol);
    set_vcm0(v,N_mol);

    double a_prev[3 * N_mol];
    compila_matr(r, a_prev, N_mol, fLJ, var_ad);//inizializzo accelerazioni come date dal potenziale
    
    
    double r_mod, v_mod;

    E = 0; K = 0;
    for (int n = 0; n < N; ++n) {//passi temporali
        E = E*(n+1);//energia totale @t
        K = K*(n+1);//energia cinetica @t
        T = 0;//temp @t
        for (int j = 0; j < N_mol; ++j) {//particelle
            t = t0 + n * h;
            if (vel_verlet(t, r, v, h, j, 3, a_prev, fLJ, var_ad)) {printf("ERRORE");}

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
double fLJ(double *r, double *args, int indice){//arg[0]=eps, arg[1]=sigma
    if((r[indice])>=args[2]/2){//manda la forza a zero se piu distante di L/2
        return 0;
    }
    else{
        double r6=r[indice]*r[indice]*r[indice]*r[indice]*r[indice]*r[indice];
        double sigma6=args[1]*args[1]*args[1]*args[1]*args[1]*args[1];
        double f=4*args[0]*(12*sigma6*sigma6/(r6*r6*r[indice])+6*sigma6/(r6*r[indice]));
        return f;
    }
}

