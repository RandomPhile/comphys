#include <iostream>
#include <cmath>
#include <fstream>
#include "eq_diff.h"

using namespace std;

ofstream dati;

double f(double r, double P, double m, double *fargs);
double g(double r, double P, double m, double *gargs);

int main() {
    dati.open("dati.dat");
    int N = 1e4;
    double Pc_max = pow(20,1), Pc_min = 1e-3;
    double r, r0 = 1e-5, m, rho;
    double M[N][3], R[N][3];
    double K[] = {0.05, 0.1, 0.01};
    double gam[] = {5.0 / 3, 4.0 / 3, 2.54};
    double cost[] ={0.414, 0.213, 0.569};


    double h = 1e-6*N;
    double args[2];
    double P;
    int j = 0;
    /*for (double Pc = Pc_min; Pc < Pc_max; Pc+=(Pc_max-Pc_min)/N) {
        P = Pc;
        m = 0;
        args[0] = gam[0]; args[1] = K[0];
        r = r0;
        while (P > 0) {
            if(j>0){
                M[j-1][0] = m;
                R[j-1][0] = r;
            }
            r += h;
            
            if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
            rho = pow(P / (args[1] * (args[0] - 1)), 1.0 / args[0]);
            //cout << r << "\t" << P << "\t" << m << "\t" << rho << endl;
            if(j==1){
                dati << r << "\t" << P << "\t" << m << "\t" << rho << endl;
                cout << r << "\t" << P << "\t" << m << "\t" << rho << endl;
            }
        }
        if(j==1){
            dati << "\n\n";
            cout<< "\n\n";
        }
        //printf("R = %f\n", r);
        
        j++;
    }
    
    
    j = 0;
    for (double Pc = Pc_min; Pc < Pc_max; Pc+=(Pc_max-Pc_min)/N) {
        P = Pc;
        m = 0;
        args[0] = gam[1]; args[1] = K[1];
        r = r0;
        while (P > 0) {
            if(j>0){
                M[j-1][1] = m;
                R[j-1][1] = r;
            }
            r += h;
            
            if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
            rho = pow(P / (args[1] * (args[0] - 1)), 1.0 / args[0]);
            //cout << r << "\t" << P << "\t" << m << "\t" << rho << endl;
            if(j==1){
                dati << r << "\t" << P << "\t" << m << "\t" << rho << endl;
                cout << r << "\t" << P << "\t" << m << "\t" << rho << endl;
            }
        }
        if(j==1){
            dati << "\n\n";
            cout<< "\n\n";
        }
        //printf("R = %f\n", r);
        j++;
    }
    
    for(int t=0; t<2; t++){
        for (int j = 0; j < N; j++) {
            if(R[j][t]<1.5 && R[j][t]>0.15){
                if(fabs(pow(M[j][t],2-gam[t])*pow(R[j][t],3*gam[t]-4)-cost[t])<h){
                    dati << M[j][t] << "\t" << R[j][t] << "\t" << pow(M[j][t],2-gam[t])*pow(R[j][t],3*gam[t]-4) << endl;
                    cout << M[j][t] << "\t" << R[j][t] << "\t" << pow(M[j][t],2-gam[t])*pow(R[j][t],3*gam[t]-4) << endl;
                }
            }
        }
        if(t<1){
            dati<<"\n\n";
            cout<<"\n\n";
        }
    }*/
    Pc_min=1e-5;
    Pc_max=1500;
    
    for (double Pc = Pc_min; Pc < Pc_max; Pc+=(Pc_max-Pc_min)/N) {
        P = Pc;
        m = 0;
        args[0] = gam[2]; args[1] = K[2];
        r = r0;
        while (P > 0) {
            if(j>0){
                M[j-1][2] = m;
                R[j-1][2] = r;
            }
            r += h;
            
            if (runge_kutta(r, &P, &m, h, f, args, g, args)) {printf("ERRORE");}
            rho = pow(P / (args[1] * (args[0] - 1)), 1.0 / args[0]);
            //cout << r << "\t" << P << "\t" << m << "\t" << rho << endl;
            if(j==1){
                dati << r << "\t" << P << "\t" << m << "\t" << rho << endl;
                //cout << r << "\t" << P << "\t" << m << "\t" << rho << endl;
            }
        }
        if(j==1){
            dati << "\n\n";
            //cout<< "\n\n";
        }
        //printf("R = %f\n", r);
        
        j++;
    }
    for(int t=2; t<3; t++){
        for (int j = 0; j < N; j++) {
            if(R[j][t]<1.5 && R[j][t]>0.15){
                if(fabs(pow(M[j][t],2-gam[t])*pow(R[j][t],3*gam[t]-4)-cost[t])<h){
                    dati << M[j][t] << "\t" << R[j][t] << "\t" << pow(M[j][t],2-gam[t])*pow(R[j][t],3*gam[t]-4) << endl;
                    cout << M[j][t] << "\t" << R[j][t] << "\t" << pow(M[j][t],2-gam[t])*pow(R[j][t],3*gam[t]-4) << endl;
                }
            }
        }
        if(t<1){
            dati<<"\n\n";
            cout<<"\n\n";
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
