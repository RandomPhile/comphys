#include <iostream>
#include <cmath>
#include <fstream>

#include "header.h"
//#include "funzioni.h"
#include "es11.h"
#include "distribuzioni.h"

/* NOTE:
-se K o V esplodono molto probabilmente dt è troppo grande
-per stimare sigma t.c. T=1.1, uso N grande ()
-M=1, caso=2:  non esiste un sigma per far convergere T a T_req=1.1
*/

/*** variabili globali ***/
//CC, BCC, FCC
int M = 2; //1,2,4
int N = M * pow(4, 3); //numero di particelle
ofstream dati, gnuplot, risultati;

struct coppia {
    double rho;
    double sigma;
};






int main() {
    srand(1);//default seed = 1
    double dt      = 0.01;//passo temporale
    double t1      = 20;//durata simulazione

    coppia coppie[] = {//aggiornate con M=2,n=6,dt=0.01,t1=20 (6+6 minuti)
        {.rho = 0.01, .sigma = 1.04332},
        {.rho = 0.1, .sigma = 0.828379},
        {.rho = 0.2, .sigma = 0.670309},
        {.rho = 0.3, .sigma = 0.656069},
        {.rho = 0.4, .sigma = 0.762703},
        {.rho = 0.5, .sigma = 0.925801},
        {.rho = 0.6, .sigma = 1.12905},
        {.rho = 0.7, .sigma = 1.28666},
        {.rho = 0.8, .sigma = 1.41451},
        {.rho = 0.9, .sigma = 1.45276},
        {.rho = 1.0, .sigma = 1.38634},
        {.rho = 1.1, .sigma = 1.39870},
        {.rho = 1.15, .sigma = 1.41261},
        {.rho = 1.2, .sigma = 1.45959}
    };

    int caso_min = 0;//mettere 0 per avere P(rho)
    int caso_max;

    if (caso_min == 0) {
        caso_max = sizeof(coppie) / sizeof(*coppie);
    } else {
        caso_max = caso_min + 1;
    }
    //###################################################
    struct vec r[N];
    int reticolo   = log2(M);
    double T_req = 1.1;

    risultati.open("risultati.dat");
    gnuplot.open("gnuplot.dat");
    gnuplot << caso_min << endl;
    gnuplot.close();

    dati.open("dati.dat");
    for (int caso = caso_min; caso < caso_max; ++caso) {
        double L = cbrt(N / coppie[caso].rho);
        double r_c = L / 2;

        double Delta=L/100;

        double t = 0;
        int N_t = (t1 - t) / dt;

        crea_reticolo(r, L);


        double E = 0, T = 0, P = 0;
        double K = 0, V = 0, W = 0;
        double K_c = 0, V_c = 0, W_c = 0;

        for (int i = 0; i < N_t; ++i) {//tempo
            V_c = V_c * (i + 1.0);
            W_c = W_c * (i + 1.0);
            
            MRT2(r, &V, &W, N, Delta, T_req, L, r_c);
            
            V_c = (V_c + V) / (i + 2.0);
            W_c = (W_c + W) / (i + 2.0);

            P = (1 + W_c / (3.0 * T_req)); //P su rho*k_B*T_req
            // P = coppie[caso].rho * (1 + W_c / (3.0 * T_req)); //P su k_B*T_req

            if (caso_min != 0) {
                dati << t << "\t" << V << "\t" << P << endl;
            }
            t += dt;
        }
        LOG(P)
        if (caso_min == 0) {
            dati << coppie[caso].rho << "\t" << P << endl;
        }
        cout << "Rho = " << coppie[caso].rho << "\t\tT = " << T_req << "\nSigma = " << coppie[caso].sigma << "\t->\t" << coppie[caso].sigma*sqrt(T_req / T) << "\n" << endl;
        risultati << "Rho = " << coppie[caso].rho << "\t\tT = " << T_req << "\nSigma = " << coppie[caso].sigma << "\t->\t" << coppie[caso].sigma*sqrt(T_req / T) << "\n" << endl;
    }
    dati.close();
    risultati.close();

    return 0;
}

