#include <iostream>
#include <cmath>
#include <fstream>
#include "eq_diff.h"
ofstream dati, dati1;

double f(double r, double P, double m, double *fargs);
double g(double r, double P, double m, double *gargs);

int main() {
    dati.open("dati.dat"); dati1.open("dati1.dat");

    double K[] = {0.05, 0.1, 0.01};
    double gam[] = {5.0 / 3, 4.0 / 3, 2.54};

    double m, P, r, rho, args[2], cost;
    double r0 = 0.0001, h = 1e-2;
    double Pc_min = 5, Pc_max = 30;

    int N = 10;//da scrivere anche su plot.plt

<<<<<<< HEAD
    int N = 1000;//da scrivere anche su plot.plt
=======
>>>>>>> parent of e18495e (Trovato range di pressioni Pc)
    double R[3][N], M[3][N], Pc[3][N];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < N; ++j) {
            dati << "\n\n(i,j) = (" << i << "," << j << ")" << endl;
            //dati << "r\tP\tm\trho" << endl;
            Pc[i][j] = Pc_min + (Pc_max - Pc_min) * j;
            P = Pc[i][j];
            m = 0;
            r = r0;

            args[0] = gam[i]; args[1] = K[i];

            while (P > 0) {
                rho = pow(P / (args[1] * (args[0] - 1)), 1.0 / args[0]);
                dati << r << "\t" << P << "\t" << m << "\t" << rho << endl;
                R[i][j] = r; M[i][j] = m;
                runge_kutta(r, &P, &m, h, f, args, g, args);
                r += h;
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        dati1 << "\n\nStella " << i << endl;
        for (int j = 0; j < N; ++j) {
<<<<<<< HEAD
            if(R[i][j]*R0<=50 && R[i][j]*R0>=3){
                cost = pow(M[i][j], 2 - gam[i]) * pow(R[i][j], 3 * gam[i] - 4);
                dati1 << Pc[i][j] << "\t" << R[i][j]*R0 << "\t" << M[i][j] << "\t" << cost << "\t" << endl;
            }
=======
            cost = pow(M[i][j], 2 - gam[i]) * pow(R[i][j], 3 * gam[i] - 4);
            dati1 << Pc[i][j] << "\t" << R[i][j] << "\t" << M[i][j] << "\t" << cost << "\t" << endl;
>>>>>>> parent of e18495e (Trovato range di pressioni Pc)
        }
    }

    dati.close(); dati1.close();
    return 0;
}

double f(double r, double P, double m, double *fargs) {
    return - (m * pow(P / (fargs[1] * (fargs[0] - 1)), 1 / fargs[0]) / pow(r, 2)) ;
}

double g(double r, double P, double m, double *gargs) {
    return pow(r, 2) * pow(P / (gargs[1] * (gargs[0] - 1)), 1 / gargs[0]);
}
