#include "funzioni.h"

const bool relativ = 0;//0 o 1 per la correzione relativistica
const bool errore = 1;//0 o 1 per la parte dell'errore (solo non relativistica)
const int N = 20;//numero valori Pc

int main() {
    /*** costanti simulazione ***/
    const double h = 0.001, r0 = 1e-5;

    double Pc_min[] = {6e-6, 3e-2, 4e-7};
    double Pc_max[] = {1e7, 3e3, 2e5};

    /*** costanti problema ***/
    double K[]   = {0.05, 0.1, 0.01};
    double gam[] = {5.0 / 3, 4.0 / 3, 2.54};
    const double R0    = 20.06145;//km
    const double M0    = 13.6557594;//ums

    /*** variabili ***/
    double **R, **M, **Pc;
    R = new double*[3]; M = new double*[3]; Pc = new double*[3];
    for (int i = 0; i < 3; i++) {
        R[i] = new double[N];
        M[i] = new double[N];
        Pc[i] = new double[N];
    }

    if (errore == 0) {
        risolvi_stelle(r0, h, gam, K, Pc_max, Pc_min, M, Pc, R);
        if (relativ) {
            stampa_valori_rel(M, Pc, R, R0, M0, gam);
        } else {
            stampa_valori(M, Pc, R, R0, M0, gam);
        }
    } else {
        
        /*** parte dell'errore ***/

        int N_E = 12; //numero punti h_E
        double Pc_E = 47.367; //una pressione esempio
        double h_E = 0.001; //step di partenza
        double args[] = {gam[0], K[0]};

        double m1, P1, r, m2, P2;
        double epsilon1, epsilon2, P_old1 = Pc_E, P_old2 = Pc_E;
        int count;

        ofstream errore; errore.open("out/errore.dat");
        for (int i = 0; i < N_E; ++i) {

            /*inizializzo variabili*/
            epsilon1 = 0; epsilon2 = 0;
            P1 = Pc_E; m1 = 0; r = r0;
            P2 = Pc_E; m2 = 0;
            count = 0;
            while (P1 > 0 &&  P2 > 0) {
                if (r != r0 && P1 > 0 && P2 > 0) {
                    count++;
                    epsilon1 += fabs((P1 - P_old1) / P_old1);
                    epsilon2 += fabs((P2 - P_old2) / P_old2);
                }

                P_old1 = P1; P_old2 = P2;
                runge_kutta(r, &P1, &m1, h_E, dP, args, dm, args);
                eulero_exp(r, &P2, &m2, h_E, dP, args, dm, args);

                r += h_E;
            }
            epsilon1 /= count;//per la media di epsilon
            epsilon2 /= count;
            cout << h_E << "\t\t" << epsilon1 << "\t" << epsilon2 << endl;
            errore << h_E << "\t" << epsilon1 << "\t" << epsilon2 << endl;

            h_E /= 2;
        }
        errore.close();
    }
    plot();
    return 0;
}
