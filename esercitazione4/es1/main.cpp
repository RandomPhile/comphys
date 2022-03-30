#include <iostream>
#include <cmath>
#include <fstream>
#include "eq_diff.h"
ofstream dati;

double f(double r, double P, double m, double *fargs);
double g(double r, double P, double m, double *gargs);

int main() {
    dati.open("dati.dat");
    int N=10;
    double r=0,r0 = 1e-10, h, m = 1e-6, rho, P0=5;
    double M[N][3], R[N][3];
    h = 1e-5;
    for(int j=0; j<N; j++){
        P0=j*4+P0;
        double P = P0;

        double gam= 5.0 / 3, K = 0.05;; //gas fermi non rel
        double args[] = {gam, K};

        
        r0 = 1e-5;
        P=P0;
        m=0;
        
        int i = 0;
        while (P > 0) {
            if(m!=NAN)
                M[j][0]=m;
            r = r0 + i * h;
            if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
            rho=pow(P / (args[1] * (args[0] - 1)), 1 / args[0]);
            dati << r << "\t" << P << "\t" << m << "\t"<<rho<<endl;
            i++;
        }
        R[j][0]=r;
        //printf("R = %f\n", r);
        dati<<"\n\n";
        
        r0 = 1e-5;
        P=P0;
        m=0;
        
        gam=4.0/3;
        K=0.1;
        args[0]=gam;
        args[1]=K;
        i = 0;
        while (P > 0) {
            if(m!=NAN)
                M[j][1]=m;
            r = r0 + i * h;
            if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
            rho=pow(P / (args[1] * (args[0] - 1)), 1 / args[0]);
            dati << r << "\t" << P << "\t" << m << "\t"<<rho<<endl;
            i++;
        }
        R[j][1]=r;
        //printf("R = %f\n", r);
        dati<<"\n\n";
        
        
        r0 = 1e-5;
        P=P0;
        m=0;
        
        gam=2.54;
        K=0.01;
        args[0]=gam;
        args[1]=K;
        i = 0;
        while (P > 0) {
            if(m!=NAN)
                M[j][2]=m;
            r = r0 + i * h;
            if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
            rho=pow(P / (args[1] * (args[0] - 1)), 1 / args[0]);
            dati << r << "\t" << P << "\t" << m << "\t"<<rho<<endl;
            i++;
        }
        R[j][2]=r;
        //printf("R = %f\n", r);
        dati<<"\n\n";
    
    }
    for(int t=0; t<3; t++){
        for(int j=0; j<N; j++){
            dati<<M[j][t]<<"\t"<<R[j][t]<<endl;
            cout<<M[j][t]<<"\t"<<R[j][t]<<endl;
        }
        if(t<2){
            cout<<"\n\n";
            dati<<"\n\n";
        }
    }
    dati.close();
    return 0;
}

double f(double r, double P, double m, double *fargs) {
    return - (m * pow(P / (fargs[1] * (fargs[0] - 1)), 1 / fargs[0])/ pow(r, 2)) ;
}

double g(double r, double P, double m, double *gargs) {
    return pow(r, 2) * pow(P / (gargs[1] * (gargs[0] - 1)), 1 / gargs[0]);
}
