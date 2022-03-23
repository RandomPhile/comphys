#include <iostream>
#include <fstream>
#include <cmath>
#include "eq_diff.h"

using namespace std;
ofstream dati;

//putatore NULL
double f(double t, double x, double y, double *fargs);
double g(double t, double x, double y, double *gargs);

double lam;

int main() {
    dati.open("dati.dat");
    double t ,t0 = 0, h, x=1, y=0,E,T;
    T = 1;
    double t1 = 25*T;
    lam = 1/T;
    int N = 10000;
    h = (t1-t0)/( (double) N);
    //printf("t\t\t\tx\t\t\ty\n");
    
    for (int n = 0; n < N; n++) {
        t = t0 + n*h;
        if (eulero_exp(t, &x, &y, h, f, NULL, g, NULL)){printf("ERRORE");}
        E=0.5*pow(x,2)+0.5*pow(y,2);
        //printf("%f\t%f\t%f\t%f\n", t,x,y,E);
        dati<<t<<"\t"<<x<<"\t"<<y<<"\t"<<E<<endl;
    }
    dati<<"\n\n\n";
    //riporto in condizione iniziale
    t0 = 0;
    t1 = 25;
    x=1;
    y=0;
    for (int n = 0; n < N; n++) {
        t = t0 + n*h;
        if (eulero_cromer(t, &x, &y, h, f, NULL, g, NULL)){printf("ERRORE");}
        E=0.5*pow(x,2)+0.5*pow(y,2);
        //printf("%f\t%f\t%f\t%f\n", t,x,y,E);
        dati<<t<<"\t"<<x<<"\t"<<y<<"\t"<<E<<endl;
    }
    
    dati<<"\n\n\n";
    //riporto in condizione iniziale
    t0 = 0;
    t1 = 25;
    x=1;
    y=0;
    for (int n = 0; n < N; n++) {
        t = t0 + n*h;
        if (eulero_imp(t, &x, &y, h, f, NULL, g, NULL)){printf("ERRORE");}
        E=0.5*pow(x,2)+0.5*pow(y,2);
        //printf("%f\t%f\t%f\t%f\n", t,x,y,E);
        dati<<t<<"\t"<<x<<"\t"<<y<<"\t"<<E<<endl;
    }

    dati<<"\n\n\n";
    //riporto in condizione iniziale
    t0 = 0;
    t1 = 25;
    x=1;
    y=0;
    for (int n = 0; n < N; n++) {
        t = t0 + n*h;
        if (runge_kutta(t, &x, &y, h, f, NULL, g, NULL)){printf("ERRORE");}
        E=0.5*pow(x,2)+0.5*pow(y,2);
        //printf("%f\t%f\t%f\t%f\n", t,x,y,E);
        dati<<t<<"\t"<<x<<"\t"<<y<<"\t"<<E<<endl;
    }
    dati.close();
    return 0;
}

double f(double t, double x, double y, double *fargs){
    return y;
};
double g(double t, double x, double y, double *gargs){
    return -x;
};

