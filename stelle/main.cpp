#include <iostream>
#include <cmath>
#include <fstream>
#include "eq_diff.h"

//variabili globali
ofstream dati, rect, gnuplot;
bool relativ = 1;//0 o 1 per la correzione relativistica
bool errore = 0;//0 o 1 per la parte dell'errore (solo non relativistica)
int N = 20;

//giocare col valore di h (0.001) e modificare la funzione runge-kutta o eulero_exp
int main() {
    double h = 0.0001, r0 = 1e-5;
    bool funzione = 0;//0: runge-kutta, 1: eulero_exp

    double K[] = {0.05, 0.1, 0.01};
    double gam[] = {5.0 / 3, 4.0 / 3, 2.54};
    double R0 = 20.06145;//km
    double M0 = 13.6557594;//ums
    double Pc_min[] = {6e-6, 3e-2, 4e-7};
    double Pc_max[] = {1e7, 3e3, 2e5};

    double **R, **M, **Pc;
    R = new double*[3]; M = new double*[3]; Pc = new double*[3];
    for (int i = 0; i < 3; i++) {
        R[i] = new double[N];
        M[i] = new double[N];
        Pc[i] = new double[N];
    }

    dati.open("dati.dat");
    risolvi_stelle(r0, h, gam, K, Pc_max, Pc_min, M, Pc, R, funzione);
    if (relativ) {
        stampa_valori_rel(M, Pc, R, R0, M0, gam);
    } else {
        stampa_valori(M, Pc, R, R0, M0, gam);
    }
    if (errore) {
        risolvi_stelle(r0, h / 2, gam, K, Pc_max, Pc_min, M, Pc, R, funzione);
        stampa_valori(M, Pc, R, R0, M0, gam);
    }
    dati.close();

    gnuplot.open("gnuplot.dat");
    gnuplot << relativ << "\n" << errore << endl;
    gnuplot.close();
    return 0;
}
