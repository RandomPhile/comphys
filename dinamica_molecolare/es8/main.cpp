#include "header.h"
#include "funzioni.h"
#include "es7.h"

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
    double t1      = 10;//durata simulazione

    coppia coppie[] = {
        {.rho = 0.01,.sigma = 1.04435},
        {.rho = 0.1, .sigma = 0.84084},
        {.rho = 0.2, .sigma = 0.724991},
        {.rho = 0.3, .sigma = 0.715884},
        {.rho = 0.4, .sigma = 0.781342},
        {.rho = 0.5, .sigma = 0.871225},
        {.rho = 0.6, .sigma = 1.08453},
        {.rho = 0.7, .sigma = 1.21159},
        {.rho = 0.8, .sigma = 1.37527},
        {.rho = 0.9, .sigma = 1.4332},
        {.rho = 1.0, .sigma = 1.49973},
        {.rho = 1.1, .sigma = 1.04971},
        {.rho = 1.2, .sigma = 1.44088}
    };

    // double rho[]   = {0.01, 0.1, 0.8, 1.2};
    // double sigma[] = {1.04, 0.863, 1.4, 1.4558}; //BCC

    int caso_min = 0;//mettere 0 per avere P(rho)
    int caso_max;

    if (caso_min == 0) {
        caso_max = sizeof(coppie) / sizeof(*coppie);
    } else {
        caso_max = caso_min+1;
    }
    //###################################################
    struct vec r[N], v[N], a[N];
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

        double t = 0;
        int N_t = (t1 - t) / dt;

        crea_reticolo(r, L);
        crea_reticolo(a, 0);
        distr_gauss(v, coppie[caso].sigma);

        vec *dr[N]; vec_2D(dr, N);
        a_LJ(r, a, dr, r_c, L);


        double E = 0, T = 0, P = 0;
        double K = 0, V = 0, W = 0;
        double K_c = 0, V_c = 0, W_c = 0;

        for (int i = 0; i < N_t; ++i) {//tempo
            K_c = K_c * (i + 1.0);
            V_c = V_c * (i + 1.0);
            W_c = W_c * (i + 1.0);

            vel_verlet(r, v, a, dt, r_c, L, &K, &V, &W);
            v_cm_0(v);
            K_c = (K_c + K) / (i + 2.0);
            V_c = (V_c + V) / (i + 2.0);
            W_c = (W_c + W) / (i + 2.0);

            T = 2.0 * K_c / (3.0 * N);
            P = coppie[caso].rho * (1 + W_c / (3.0 * T_req)); //P su T_req

            E = K + V;
            if (caso_min != 0) {
                dati << t << "\t" << K << "\t" << V << "\t" << E << "\t" << T << "\t" << P << endl;
            }
            t += dt;
        }
        if (caso_min == 0) {
            dati << coppie[caso].rho << "\t" << P << endl;
        }
        cout << "Rho = " << coppie[caso].rho << "\t\tT = " << T << "\nSigma = " << coppie[caso].sigma << "\t->\t" << coppie[caso].sigma*sqrt(T_req / T) << "\n" << endl;
        risultati << "Rho = " << coppie[caso].rho << "\t\tT = " << T << "\nSigma = " << coppie[caso].sigma << "\t->\t" << coppie[caso].sigma*sqrt(T_req / T) << "\n" << endl;
    }
    dati.close();
    risultati.close();

    return 0;
}

