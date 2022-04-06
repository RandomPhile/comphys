#include <iostream>
#include <fstream>
#include <cmath>
#include "eq_diff.h"

using namespace std;
ofstream dati;

//putatore NULL
double f(double t, double x, double y, double *fargs);
double g(double t, double x, double y, double *gargs);
double a(double *r, int indice, double *args);

double lam;

int main() {
    dati.open("dati.dat");
    double t , t0 = 0, h, x = 1, y = 0, E, T;
    T = 1;
    double t1 = 25 * T;
    lam = 1 / T;
    int N = 10000;
    h = (t1 - t0) / ( (double) N);
    //printf("t\t\t\tx\t\t\ty\n");

    //riporto in condizione iniziale
    t0 = 0;
    t1 = 25 * T;
    double r[] = {1,1,1};
    double r_mod = sqrt(pow(r[0],2)+pow(r[1],2)+pow(r[2],2));
    double v[] = {0,0,0};
    double v_mod = sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
    double a_prev[3];
    for (int i = 0; i < 3; ++i) {
        a_prev[i] = r[i];
    }
    for (int n = 0; n < N; n++) {
        t = t0 + n * h;

        if (vel_verlet(t, r, v, h, a_prev, a, NULL)) {printf("ERRORE");}
        
        r_mod = sqrt(pow(r[0],2)+pow(r[1],2)+pow(r[2],2));
        v_mod = sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));

        E = 0.5 * pow(r_mod, 2) + 0.5 * pow(v_mod, 2);
        dati << t << "\t" << r[0] << "\t" << v[0] << "\t" << E << endl;
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

