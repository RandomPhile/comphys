#include <iostream>
#include <cmath>
#include <fstream>
#include "eq_diff.h"
ofstream dati;

double f(double r, double P, double m, double *fargs);
double g(double r, double P, double m, double *gargs);

int main() {
    dati.open("dati.dat");
    int N = 10; 
    double Pc_max = pow(50,3), Pc_min = 1e-3;
    double r, r0 = 1e-5, m, rho, Pc = 5;
    double M[N], R[N];
    double K[] = {0.05, 0.1, 0.01};
    double gam[] = {5.0 / 3, 4.0 / 3, 2.54};


    double h = 1e-3;
    double args[2];
    double P;
    int j = 0;
    for (double Pc = Pc_min; Pc < Pc_max; Pc+=(Pc_max-Pc_min)/N) {
        P = Pc;
        m = 0;
        args[0] = gam[0]; args[1] = K[0];
        r = r0;
        while (P > 0) {
            M[j] = m;
            r += h;
            
            if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
            rho = pow(P / (args[1] * (args[0] - 1)), 1.0 / args[0]);
            //cout << r << "\t" << P << "\t" << m << "\t" << rho << endl;
            dati << r << "\t" << P << "\t" << m << "\t" << rho << endl;
        }
        R[j] = r;
        //printf("R = %f\n", r);
        dati << "\n\n";
        j++;
    }

    for (int j = 0; j < N; j++) {
        dati << M[j] << "\t" << R[j] << "\t" << pow(M[j],2-gam[0])*pow(R[j],3*gam[0]-4) << endl;
        cout << M[j] << "\t" << R[j] << "\t" << pow(M[j],2-gam[0])*pow(R[j],3*gam[0]-4) << endl;
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
