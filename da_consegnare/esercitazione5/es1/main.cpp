#include <iostream>
#include <cmath>
#include <fstream>
#include "eq_diff.h"
ofstream dati, dati1;

double dP(double r, double P, double m, double *args);
double dm(double r, double P, double m, double *args);

int main() {
    dati.open("dati.dat"); dati1.open("dati1.dat");

    double K[] = {0.05, 0.1, 0.01};
    double gam[] = {5.0 / 3, 4.0 / 3, 2.54};

    double m, P, r, args[4], cost;
    double r0 = 0.00001, h = 1e-3, h_pc;
    
    // double Pc_min[] = {6e-6, 3e-2, 4e-7};
    // double Pc_max[] = {1e5, 3e3, 2e5};//ho cambiato per stella1 da 1e7 a 1e5

    //le cose interessanti succedono tra:
    double Pc_min[] = {6e-6, 3e-2, 4e-7};
    double Pc_max[] = {1e1, 1e1, 1e1};

    int N = 200;//da scrivere anche su plot.plt
    double derivata[3][N-1];
    double R[3][N], M[3][N], Pc[3][N];

    double R0 = 20.06145;//km
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < N; ++j) {
            dati << "\n\n(i,j) = (" << i << "," << j << ")" << endl;
            //dati << "r\tP\tm\trho" << endl;
            h_pc = ((Pc_max[i] - Pc_min[i]) / (N - 1)) * j;
            Pc[i][j] = Pc_min[i] + h_pc;
            P = Pc[i][j];//P iniziale
            m = 0;
            r = r0;

            args[0] = gam[i]; args[1] = K[i];

            while (P > 0) {
                args[2] = pow(P / (args[1] * (args[0] - 1)), 1.0 / args[0]);//rho
                args[3] = args[2] + K[i] * pow(args[2], gam[i]);//epsilon
                dati << r << "\t" << P << "\t" << m << "\t" << args[2] << endl;
                R[i][j] = r; M[i][j] = m;
                runge_kutta(r, &P, &m, h, dP, args, dm, args);
                r += h;
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        dati1 << "\n\nStella " << i << endl;
        for (int j = 0; j < N; ++j) {
            cost = pow(M[i][j], 2 - gam[i]) * pow(R[i][j], 3 * gam[i] - 4);
            dati1 << Pc[i][j] << "\t" << R[i][j]*R0 << "\t" << M[i][j] << "\t" << cost << "\t" << endl;
        }
        for (int j = 0; j < N-1; ++j) {
            derivata[i][j] = -fabs((M[i][j+1]-M[i][j])/(Pc[i][j+1]-Pc[i][j]));
        }
        int indice = distance(derivata[i],max_element(derivata[i], derivata[i] + N-1));
        cout<<Pc[i][indice]<<endl;
    }
    

    dati.close(); dati1.close();
    return 0;
}

double dP(double r, double P, double m, double *args) {
    return -((P + args[3]) * (m + pow(r, 3) * P)) / (pow(r, 2) - 2 * m * r);
}

double dm(double r, double P, double m, double *args) {
    return pow(r, 2) * args[3];
}
