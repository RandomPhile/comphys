#include <iostream>
#include <cmath>
#include <fstream>
#include "eq_diff.h"
ofstream dati;

double f(double r, double P, double m, double *fargs);
double g(double r, double P, double m, double *gargs);

int main() {
    dati.open("dati.dat");
    int N = 5;
    double r, r0 = 1e-5, m, rho, Pc[3] = {10,10,10};
    double M[N][3], R[N][3];
    double K[] = {0.05, 0.1, 0.01};
    double gam[] = {5.0 / 3, 4.0 / 3, 2.54};


    double h = 1e-5;
    double args[2];
    double P;
    int i;
    for (int j = 0; j < N; j++) {
        for (int t = 0; t < 3; ++t) {
            r = r0; m = 0;
            Pc[t] += j*4;
            P = Pc[t];
            args[0] = gam[t]; args[1] = K[t];

            i = 0;
            while (P > 0) {
                if (m != NAN) {
                    M[j][t] = m;
                    r += i * h;
                } else {break;}
                if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
                rho = pow(P / (args[1] * (args[0] - 1)), 1.0 / args[0]);
                dati << r << "\t" << P << "\t" << m << "\t" << rho << endl;
                i++;
                
            }
            R[j][t] = r;
            //printf("R = %f\n", r);
            dati << "\n\n";
        }
    }

    for (int t = 0; t < 3; t++) {
        for (int j = 0; j < N; j++) {
            dati << M[j][t] << "\t" << R[j][t] << endl;
        }
        if (t < 2) {
            dati << "\n\n";
        }
    }

    dati.close();
    return 0;
}

double f(double r, double P, double m, double *fargs) {
    return - (m * pow(P / (fargs[1] * (fargs[0] - 1)), 1 / fargs[0]) / pow(r, 2)) ;
}

double g(double r, double P, double m, double *gargs) {
    return pow(r, 2) * pow(P / (gargs[1] * (gargs[0] - 1)), 1 / gargs[0]);
}
