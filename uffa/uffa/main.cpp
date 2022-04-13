#include <iostream>
#include <fstream>
#include <cmath>
//#include "eq_diff.h"

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
    int N_P = 40/0.1, N_M = 100000; //N_P PARI
    h = (t1 - t0) / ( (double) N_P);
    //printf("t\t\t\tx\t\t\ty\n");
    double r[N_P];
    //creo la distribuziuone di velocità
    if (N_M % 2 != 0) {
        N_M++;
    }
    double v[N_M];
    double a_prev[N_M];

    for (int i = 0; i < N_P; i++) {
        r[i] = 0;
    }



    for (int i = 0; i < N_M; i += 2) {
        double x1, x2;
        x1 = rand() / ((double)RAND_MAX + 1.0);
        x2 = rand() / ((double)RAND_MAX + 1.0);
        v[i] = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
        v[i + 1] = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1);
        //dati << i << "\t" << v[i] << endl;
    }
    //verifico che sia gaussiana
    for (double j = -20; j < 20; j += 40.0 / N_P) {
        int cont = 0;
        for (int i = 0; i < N_M; i++) {
            if (v[i] < j + 40.0 / N_P && v[i] > j) {
                cont++;
            }
        }
        dati << j << "\t" << cont << endl;
    }
    // }
    // for(int i=0;i<N_M;i++){
    //     for(int j=0;j<N_P;j++){
    //         if(v[i]<-(N_P/2)+(j+0.1) && v[i]>-(N_P/2)+j){
    //             r[j]++;
    //             dati<<r[j]<<"\t"<<v[i]<<endl;
    //             //cout<<r[j]<<"\t"<<v[i]<<endl;
    //             break;
    //         }
    //     }
    // }
    // double v_cm=0;
    // for(int i=0;i<N_M;i++){
    //     v_cm+=v[i];
    // }
    // v_cm/=N_M;
    // cout<<v_cm<<endl;
    // for(int i=0;i<N_M;i++){
    //     v[i]-=v_cm;
    // }
    // v_cm=0;
    // for(int i=0;i<N_M;i++){
    //     v_cm+=v[i];
    // }
    // v_cm/=N_M;
    // cout<<v_cm<<endl;

    /*double r_mod;
    double v_mod;

    for (int i = 0; i < 3; ++i) {
        a_prev[i] = r[i];
    }
    for (int n = 0; n < N_P; n++) {
        t = t0 + n * h;

        if (vel_verlet(t, r, v, h, a_prev, a, NULL)) {printf("ERRORE");}

        r_mod = sqrt(pow(r[0],2)+pow(r[1],2)+pow(r[2],2));
        v_mod = sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));

        E = 0.5 * pow(r_mod, 2) + 0.5 * pow(v_mod, 2);
        dati << t << "\t" << r[0] << "\t" << v[0] << "\t" << E << endl;
    }
    */
    dati.close();
    return 0;
}

double f(double t, double x, double y, double * fargs) {
    return y;
}
double g(double t, double x, double y, double * gargs) {
    return -x;
}
double a(double * r, int indice, double * args) {
    return -r[indice];
}

