#include "funzioni.h"

const bool relativ = 0;//0 o 1 per la correzione relativistica
const int N = 20;//numero valori di pressioni centrali

int main() {
    /*** costanti simulazione ***/
    const double h = 0.001, r0 = 1e-5;
    double Pc_min[] = {6e-6, 3e-2, 4e-7};
    double Pc_max[] = {1e7, 3e3, 2e5};

    /*** costanti problema ***/
    double K[]   = {0.05, 0.1, 0.01};
    double gam[] = {5.0 / 3, 4.0 / 3, 2.54};
    const double R0    = 20.06145;//km
    const double M0    = 13.6557594;//ums

    /*** crea variabili ***/
    double **R = new double*[3], **M = new double*[3], **Pc = new double*[3];
    for (int i = 0; i < 3; i++) {
        R[i] = new double[N]; M[i] = new double[N]; Pc[i] = new double[N];
    }

    /*** chiamata a funzioni ***/
    risolvi_stelle(r0, h, gam, K, Pc_max, Pc_min, M, Pc, R);
    if (relativ) {
        stampa_valori_rel(M, Pc, R, R0, M0, gam);
    } else {
        stampa_valori(M, Pc, R, R0, M0, gam);
    }

    /*** elimina variabili ***/
    for (int i = 0; i < 3; i++) {
        delete [] R[i]; delete [] M[i]; delete [] Pc[i];
    }
    delete [] R; delete [] M; delete [] Pc;

    /*** plot con gnuplot ***/
    plot();
    return 0;
}
