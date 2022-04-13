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

int main() {
    dati.open("dati.dat");
    double t , t0 = 0, h, E, Tau;
    double sigma=3;//uguale a sqrt(1/(m*beta))
    Tau = 1;
    double t1 = 25 * Tau;
    int N_P = 1000, N_M=50;
    h = (t1 - t0) / ( (double) N_P);
    //printf("t\t\t\tx\t\t\ty\n");
    double r[12];
    double v[N_M];
    double a_prev[N_M];
    
    for (int i=0; i<12; i++) {
        r[i]=0;
    }
    
    //creo la distribuziuone di velocità
    if(N_M%2!=0){
        N_P++;
    }
    //verifico che sia gaussiana
    for(int i=0;i<N_M;i+=2){
        double x1,x2;
        x1=rand()/((double)RAND_MAX+1.0);
        x2=rand()/((double)RAND_MAX+1.0);
        v[i]=sigma*sqrt(-2*log(1-x2))*cos(2*M_PI*x1);
        v[i+1]=sigma*sqrt(-2*log(1-x2))*sin(2*M_PI*x1);
    }
    for(int i=0;i<N_M;i+=2){
        if(v[i]<-1000){
            r[0]++;
            dati<<r[0]<<"/t"<<v[i]<<endl;
        }
        else if (v[i]>=-1000 && v[i]<-400){
            r[1]++;
            dati<<r[1]<<"/t"<<v[i]<<endl;
        }
        else if (v[i]>=-400 && v[i]<-100){
            r[2]++;
            dati<<r[2]<<"/t"<<v[i]<<endl;
        }
        else if (v[i]>=-100 && v[i]<-50){
            r[3]++;
            dati<<r[3]<<"/t"<<v[i]<<endl;
        }
        else if (v[i]>=-50 && v[i]<-15){
            r[4]++;
            dati<<r[4]<<"/t"<<v[i]<<endl;
        }
        else if (v[i]>=-15 && v[i]<0){
            r[5]++;
            dati<<r[5]<<"/t"<<v[i]<<endl;
        }
        else if (v[i]>=0 && v[i]<15){
            r[6]++;
            dati<<r[6]<<"/t"<<v[i]<<endl;
        }
        else if (v[i]>=15 && v[i]<50){
            r[7]++;
            dati<<r[7]<<"/t"<<v[i]<<endl;
        }
        else if (v[i]>=50 && v[i]<100){
            r[8]++;
            dati<<r[8]<<"/t"<<v[i]<<endl;
        }
        else if (v[i]>=100 && v[i]<400){
            r[9]++;
            dati<<r[9]<<"/t"<<v[i]<<endl;
        }
        else if (v[i]>=400 && v[i]<1000){
            r[10]++;
            dati<<r[10]<<"/t"<<v[i]<<endl;
        }
        else{
            r[11]++;
            dati<<r[11]<<"/t"<<v[i]<<endl;
        }
    }
    
    
    /*
    double r_mod;
    double v_mod = sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
    
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

double f(double t, double x, double y, double *fargs) {
    return y;
}
double g(double t, double x, double y, double *gargs) {
    return -x;
}
double a(double *r, int indice, double *args) {
    return -r[indice];
}

