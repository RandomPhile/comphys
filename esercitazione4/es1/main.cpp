#include <iostream>
#include <cmath>
#include <fstream>
#include "eq_diff.h"
ofstream dati, dati1;

double f(double r, double P, double m, double *fargs);
double g(double r, double P, double m, double *gargs);
double dm(double r, double P, double m, double *fargs);
double dP(double r, double P, double m, double *fargs);

int main() {
    dati.open("dati.dat"); dati1.open("dati1.dat");

    double K[] = {0.05, 0.1, 0.01};
    double gam[] = {5.0 / 3, 4.0 / 3, 2.54};

    double m, P, r, rho, args[2], cost, e;
    double r0 = 0.00001, h = 1e-3;
    // double Pc_min[] = {5,5,5};
    // double Pc_max[] = {10,10,10};

    double Pc_min[] = {6e-6, 3e-2, 4e-7};
    double Pc_max[] = {1e7, 3e3, 2e5};

    int N = 100;//da scrivere anche su plot.plt
    double R[3][N], M[3][N], Pc[3][N];

    double R0 = 20.06145;//km
    /*for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < N; ++j) {
            dati << "\n\n(i,j) = (" << i << "," << j << ")" << endl;
            //dati << "r\tP\tm\trho" << endl;
            Pc[i][j] = Pc_min[i] + ((Pc_max[i] - Pc_min[i]) / (N - 1)) * j;
            P = Pc[i][j];
            m = 0;
            r = r0;

            args[0] = gam[i]; args[1] = K[i];

            while (P > 0) {
                rho = pow(P / (args[1] * (args[0] - 1)), 1.0 / args[0]);
                dati << r << "\t" << P << "\t" << m << "\t" << rho << endl;
                R[i][j] = r; M[i][j] = m;
                runge_kutta(r, &P, &m, h, f, args, g, args);
                r += h;
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        dati1 << "\n\nStella " << i << endl;
        for (int j = 0; j < N; ++j) {
            if(R[i][j]*R0<=50 && R[i][j]*R0>=3){
                cost = pow(M[i][j], 2 - gam[i]) * pow(R[i][j], 3 * gam[i] - 4);
                dati1 << Pc[i][j] << "\t" << R[i][j]*R0 << "\t" << M[i][j] << "\t" << cost << "\t" << endl;
            }
        }
    }
*/
    
    
    //calcolo ora i risultati per le stelle relativistiche
    double Pc_min_rel[] = {6e-6, 3e-2, 4e-7};
    double Pc_max_rel[] = {1e7, 3e3, 2e5};

    int N_rel = 2;//da scrivere anche su plot.plt
    
    
    double R_rel[3][N], M_rel[3][N], Pc_rel[3][N];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < N; ++j) {
            dati << "\n\n(i,j) = (" << i << "," << j << ")" << endl;
            dati << "r\tP\tm\trho" << endl;
            Pc_rel[i][j] = Pc_min_rel[i] + ((Pc_max_rel[i] - Pc_min_rel[i]) / (N_rel - 1)) * j;
            P = Pc_rel[i][j];
            m = 0;
            r = r0;

            args[0] = gam[i]; args[1] = K[i];

            while (P > 0) {
                e = pow(P / (args[1] * (args[0] - 1)), 1.0 / args[0])+ P / (args[0] - 1);
                dati << r << "\t" << P << "\t" << m << "\t" << e << endl;
                R_rel[i][j] = r; M_rel[i][j] = m;
                runge_kutta(r, &P, &m, h, dP, args, dm, args);
                r += h;
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        //dati1 << "\n\nStella " << i << endl;
        cout << "\n\nStella Relativistica " << i << endl;
        for (int j = 0; j < N_rel; ++j) {
            //if(R_rel[i][j]*R0<=50 && R_rel[i][j]*R0>=3){
                cost = pow(M_rel[i][j], 2 - gam[i]) * pow(R_rel[i][j], 3 * gam[i] - 4);
                dati1 << Pc_rel[i][j] << "\t" << R_rel[i][j]*R0 << "\t" << M_rel[i][j] << "\t" << cost << "\t"<< endl;
                cout << Pc_rel[i][j] << "\t" << R_rel[i][j]*R0 << "\t" << M_rel[i][j] << "\t" << cost << "\t"<< endl;
            //}
        }
    }
    
    dati.close(); dati1.close();
    return 0;
}

double f(double r, double P, double m, double *fargs) {
    return - (m * pow(P / (fargs[1] * (fargs[0] - 1)), 1 / fargs[0]) / pow(r, 2)) ;
}
double g(double r, double P, double m, double *gargs) {
    return pow(r, 2) * pow(P / (gargs[1] * (gargs[0] - 1)), 1 / gargs[0]);
}


double dm(double r, double P, double m, double *dmargs){//assumo valgano le equazioni della volta scorsa
    double e=pow(P / (dmargs[1] * (dmargs[0] - 1)), 1 / dmargs[0]) + P / (dmargs[0] - 1);
    return pow(r, 2) * e;
}
double dP(double r, double P, double m, double *dPargs){
    double e=pow(P / (dPargs[1] * (dPargs[0] - 1)), 1 / dPargs[0]) + P / (dPargs[0] - 1);
    return -(P+e)*(m+pow(r, 3)*P)/(pow(r, 2)-2*m*r);
}
