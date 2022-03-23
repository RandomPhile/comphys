#include <iostream>
#include <cmath>
#include <fstream>
#include "eq_diff.h"
ofstream dati;

double f(double r, double P, double m, double *fargs);
double g(double r, double P, double m, double *gargs);

int main() {
    dati.open("dati.dat");
    double r,r0 = 0.01, h, m = 0, rho;
    double P = 1;

    double gam= 5.0 / 3, K = 0.05;; //gas fermi non rel
    double args[] = {gam, K};

    h = 1e-5;
    
    int i = 0;
    while (P > 0) {
        r = r0 + i * h;
        if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
        rho=pow(P / (args[1] * (args[0] - 1)), 1 / args[0]);
        dati << r << "\t" << P << "\t" << m << "\t"<<rho<<endl;
        i++;
    }
    printf("R = %f\n", r);
    dati<<"\n\n";
    
    r0 = 0.01;
    P=1;
    m=0;
    
    gam=4/3;
    K=0.1;
    args[0]=gam;
    args[1]=K;
    i = 0;
    while (P > 0) {
        r = r0 + i * h;
        if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
        rho=pow(P / (args[1] * (args[0] - 1)), 1 / args[0]);
        dati << r << "\t" << P << "\t" << m << "\t"<<rho<<endl;
        i++;
    }
    printf("R = %f\n", r);
    dati<<"\n\n";
    
    
    r0 = 0.01;
    P=1;
    m=0;
    
    gam=2.54;
    K=0.01;
    args[0]=gam;
    args[1]=K;
    i = 0;
    while (P > 0) {
        r = r0 + i * h;
        if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
        rho=pow(P / (args[1] * (args[0] - 1)), 1 / args[0]);
        dati << r << "\t" << P << "\t" << m << "\t"<<rho<<endl;
        i++;
    }
    printf("R = %f\n", r);
    dati<<"\n\n";
    
    dati.close();
    return 0;
}

double f(double r, double P, double m, double *fargs) {
    return - (m / pow(r, 2)) * pow(P / (fargs[1] * (fargs[0] - 1)), 1 / fargs[0]);
}

double g(double r, double P, double m, double *gargs) {
    return pow(r, 2) * pow(P / (gargs[1] * (gargs[0] - 1)), 1 / gargs[0]);
}
