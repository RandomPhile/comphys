#include <iostream>
#include <cmath>
#include <fstream>
#include "eq_diff.h"
ofstream dati;

double f(double r, double P, double m, double *fargs);
double g(double r, double P, double m, double *gargs);

int main() {
    dati.open("dati.dat");
    double r,r0 = 0.01, h, m = 0, E, T;
    double P = 1;

    double gam = 5.0 / 3, K = 0.05; //gas fermi non rel
    double args[] = {gam, K};


    h = 0.1;
    int i = 0;
    while (P > 0) {
        r = r0 + i * h;
        if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
        dati << r << "\t" << P << "\t" << m << endl;
        i++;
    }
    printf("R = %f\n", r);

    dati.close();
    return 0;
}

double f(double r, double P, double m, double *fargs) {
    return - (m / pow(r, 2)) * pow(P / (fargs[1] * (fargs[0] - 1)), 1 / fargs[0]);
}

double g(double r, double P, double m, double *gargs) {
    return pow(r, 2) * pow(P / (gargs[1] * (gargs[0] - 1)), 1 / gargs[0]);
}