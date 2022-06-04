#include <iostream>
#include <fstream>
#include <cmath>

#define ARMA_NO_DEBUG//non fa controlli e va piu veloce, da togliere solo se il codice funziona gia perfettamente

#include <armadillo>
#include <iomanip>
#include <algorithm>

#include "es11.h"

/*** variabili globali ***/
//CC, BCC, FCC
int M = 2; //1,2,4
int N = M * pow(6, 3); //numero di particelle

int numero_proposti=0;
int numero_accettati=0;
ofstream dati, gnuplot, risultati;

int main() {
    srand(1);//default seed = 1
    int N_t;//numero passi simulazione, vedi rigo 61 per il valore in base al tempo di equilibrazione stimato
    int N_b=20;//numero bin
    double T_req = 1.1;//temperatura adimensionale

    rowvec rho={0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
    rowvec passi_eq={400, 550, 1200, 1000, 2200, 2500, 2500, 2800, 3700, 3700, 3700, 3700, 3700, 3700, 3600, 2200};//tempi di equilibrazione con delta preso per avere circa 50%

    int caso_min = -1;//mettere -1 per avere P(rho)
    int caso_max;
    if (caso_min == -1) {
        caso_max = rho.size();
    } else {
        caso_max = caso_min + 1;
    }
    
    //###################################################
    mat r(N,3);
    cube dr(N,N,3);
    int reticolo   = log2(M);
    risultati.open("risultati.dat");
    gnuplot.open("gnuplot.dat");
    gnuplot << caso_min << endl;
    gnuplot.close();
    dati.open("dati.dat");
    int start;
    if(caso_min==-1){
        start=caso_min+1;
    }
    else{
        start=caso_min;
    }
    for (int caso = start; caso < caso_max; ++caso) {

        N_t=passi_eq(caso)+1e4;//faccio fare 10000 passi dopo l'equilibrazione per avere dei risultati carini

        double L = cbrt(N / rho(caso));
        double r_c = L / 2;
        double Delta;

        Delta=Calcola_Delta(L, rho, caso);

        crea_reticolo(r, L);//creo il reticolo iniziale

        for (int i = 0; i < N; ++i){//creo il cubo di distanze relative alla posizione iniziale
            for (int j = i + 1; j < N; ++j){
                for (int k = 0; k < 3; ++k){
                    dr(i,j,k) = r(i,k) - r(j,k);
                    dr(i,j,k) -= L * rint(dr(i,j,k)/L);//sposto in [-L/2,+L/2]
                }
            }
        }

        double P = 0, sigma_P;
        double V = 0, W = 0;
        double V_m = 0, W_m = 0, P_m = 0;
        double var_P = 0;

        for (int i = 0; i < N_t; ++i) {//passi
            if(i>passi_eq(caso)){
                V_m = V_m * (i-passi_eq(caso));
                W_m = W_m * (i-passi_eq(caso));
                var_P = var_P * (i - passi_eq(caso));
            }

            MRT2(r, &V, &W, N, Delta, T_req, L, r_c, dr);
            // cout<<V<<endl;

            if(i>passi_eq(caso)){//tempo di equilibrazione
                V_m = (V_m + V) / (i-passi_eq(caso) + 1.0);
                W_m = (W_m + W) / (i-passi_eq(caso) + 1.0);

                P = (1 + W / (3.0 * T_req));
                P_m = (1 + W_m / (3.0 * T_req)); //P su rho*k_B*T_req
                var_P = (var_P + (P - P_m) * (P - P_m)) / (i - passi_eq(caso) + 1.0);//calcolo la varianza "ordinaria"

                if (caso_min != -1) {
                    dati << i << "\t" << V_m << "\t" << P_m  << "\t" << V << "\t" << P << "\t" << sqrt(var_P) << endl;
                }
            }

            if(i==(int)(N_t/2) && caso_min!=-1){
                cout<<"Sono a metà dei cicli montecarlo, dai che ce la faccio"<<endl;
            }
        }

        cout << "Densità = "<< rho(caso) << "\t" <<"Pressione = " << P_m << endl;
        if (caso_min == -1) {
            dati << rho(caso) << "\t" << P << "\t" << sqrt(var_P) << endl;
        }
        else{//calcolo della gdr
            risultati << rho(caso) << "\t" << P_m << "\t" << sqrt(var_P) << endl;
            gdr_funz(r, L, rho(caso), N_b);
        }
        cout <<"Ho accettato il " <<(double)numero_accettati/(double)numero_proposti*100 << "%. I passi totali erano "<<N_t<<"\n\n";
    }

    dati.close();
    risultati.close();
    return 0;
}
void blocking(int N_t, int part_per_B, ifstream &dati, ifstream &risultati){
    int N_B = floor(N_t / part_per_B);

    double variabile_inutile, P_media, P_mB=0;
    double var_PB=0;
    rowvec P_m(N_B), P(N_t);

    dati.open("dati.dat");
    risultati.open("risultati.dat");

    for (int i = 0; i < N_t; ++i){
        dati >> variabile_inutile >> variabile_inutile >> variabile_inutile >> variabile_inutile >> P(i) >>variabile_inutile;
    }
    risultati >> variabile_inutile >> P_media >> variabile_inutile;
    dati.close();
    risultati.close();

    for (int i = 0; i < N_B; ++i){//calcolo le medie sui blocchi
        double somma = 0;
        for (int j = 0; j < part_per_B; ++j){
            somma += P(i * N_B + j);
        }
        P_m(i) = somma / part_per_B;
    }
    for (int i = 0; i < N_B; ++i){//calcolo la media complessiva come media sui blocchi
        P_mB += P_m(i) / N_B;
        var_PB += (P_m(i) - P_media) * (P_m(i) - P_media) / N_B;
    }
}
























