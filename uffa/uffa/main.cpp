#include <iostream>
#include <fstream>
#include <cmath>
#include "eq_diff.h"
#include "funzioni_utili.h"

using namespace std;
ofstream dati;

//putatore NULL
double f(double t, double x, double y, double *fargs);
double g(double t, double x, double y, double *gargs);
double a(double *r, int indice, double *args);

int main() {
    dati.open("dati.dat");
    double t , t0 = 0, h, E, Tau;
    double sigma = 3; //uguale a sqrt(1/(m*beta))
    Tau = 1;
    double t1 = 25 * Tau;
    
    int N_Partiz = 100, N_Mol = 100000; //N_P PARI
    h = (t1 - t0) / ( (double) N_Partiz);
    if (N_Mol % 2 != 0) {
        N_Mol++;
    }
    if (N_Partiz % 2 != 0) {
        N_Partiz++;
    }
    double r[N_Mol][3];
    double v[N_Mol][3];
    double a_prev[N_Mol][3];
    
    //creo la distribuziuone di velocità
    gauss_distr(v, sigma, N_Mol);
    
    
    //verifica che sia gaussiana
    //int N_bin=40/0.1;
    //verifica_gauss(v, N_bin, N_Mol);
    resetta_matr(r, N_Mol);
    //set v_cm=0 come media delle velocità
    set_vcm0(v, N_Mol);
    double r_mod;
    double v_mod;
    
    for(int j=0;j<3;j++){
        for (int i = 0; i < 3; ++i) {
            a_prev[i][j] = r[i][j];
        }
    }
//    for(int i=0;i<N_Mol;i++){
//        for (int n = 0; n < N_Partiz; n++) {
//            t = t0 + n * h;
//
//            if (vel_verlet(t, &r[3], &v[3], h, N_Mol, &a_prev[3], a, NULL)) {printf("ERRORE");}
//
//            r_mod = sqrt(pow(r[i][0],2)+pow(r[i][1],2)+pow(r[i][2],2));
//            v_mod = sqrt(pow(v[i][0],2)+pow(v[i][1],2)+pow(v[i][2],2));
//
//            E = 0.5 * pow(r_mod, 2) + 0.5 * pow(v_mod, 2);
//            //dati << t << "\t" << r[0] << "\t" << v[0] << "\t" << E << endl;
//        }
//    }
    
    dati.close();
    return 0;
}

double f(double t, double x, double y, double * fargs) {
    return y;
}
double g(double t, double x, double y, double * gargs) {
    return -x;
}
double a(double * r[3], int indice, double * args) {
    return -r[indice][1];
}


