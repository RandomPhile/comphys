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
int M = 4; //1,2,4 scegliere per la configurazione
int N = M * pow(6, 3); //numero di particelle

int numero_proposti = 0;
int numero_accettati = 0;
ofstream dati, gnuplot, risultati, dati_blocking;











// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------
// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------

double mod(cube &r, double riga, double colonna);
double V_LJ(double r, double L);
void MRT2(mat &r, double *V, double *W, int N, double Delta, double T, double L, cube &dr, cube &dr_n);
double dV_LJ(double r, double L);
void calcola_osservabili(rowvec rho, rowvec passi_eq, int caso_min, int N_t, int q, double T_req);
// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------
// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------








int main() {
    

// -----------------------------------------VARIABILI----------------------------VARIABILI------------------------------VARIABILI-----------------------------------------------
// -----------------------------------------VARIABILI----------------------------VARIABILI------------------------------VARIABILI-----------------------------------------------
    srand(1);//default seed = 1
    int N_t;//numero passi simulazione, vedi rigo 61 per il valore in base al tempo di equilibrazione stimato
    int N_b = 26;//numero bin
    double T_req = 1.1;//temperatura adimensionale

    rowvec rho = {0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2};//, 1.3, 1.4
    rowvec passi_eq = {1000, 2000, 3000, 4000, 6000, 9000, 9000, 9000, 9000, 1e4, 1e4, 1e4, 1e4, 1e4};//tempi di equilibrazione con delta preso per avere circa 50% , 3600, 2200
    
    int caso_min = 1;//mettere -1 per avere P(rho)
    int q = 1e5;//numero punti in piu rispetto al tempo di equilibrazione 
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    calcola_osservabili(rho, passi_eq, caso_min, N_t, q, T_req);

    // plot_pressioni_energie(caso_min);//plot pressioni ed energie

    blocking(passi_eq(caso_min)+q);//plot blocking
 
    jackknife(passi_eq(caso_min)+q);

    // gdr_funz(rho(caso_min), N_b);

    blocking_plot();

    return 0;
}

























// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------
// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------
double mod(cube &r, double riga, double colonna){//calcolo il modulo della posizione relativa delle particelle i e j
    double mod = sqrt(pow1(r(riga,colonna,0),2)+pow1(r(riga,colonna,1),2)+pow1(r(riga,colonna,2),2));
    return mod;
}
double V_LJ(double r, double L) {
    //potenziale va a zero in modo continuo sul bordo della scatola
    if (r < L / 2 && r != 0) {
        double VL_2 = 4 * (pow1(2 / L, 12) - pow1(2 / L, 6));
        return 4 * (pow1(1 / r, 12) - pow1(1 / r, 6)) - VL_2;
    } else {
        return 0;
    }
}
double dV_LJ(double r, double L) {//calcola la derivata del potenziale e la moltiplica per r
    if (r < L / 2 && r != 0) {
        return 24 * (pow1(1 / r, 6) - 2 * pow1(1 / r, 12));
    } else {
        return 0;
    }
}
void MRT2(mat &r, double *V, double *W, int N, double Delta, double T, double L, cube &dr, cube &dr_n){//N=N_mol
    rowvec r_n(3);
    
    double s = rand() / (RAND_MAX + 1.0);//eseguo prima il rand perche da problemi senno
    int n = rint((N-1) * s);//trovo la molecola che viene modificata da MTR2

    r_n(0) = r(n,0) + Delta * ((rand() / (RAND_MAX + 1.0)) - 0.5);//modifico le posizioni della particella n
    r_n(1) = r(n,1) + Delta * ((rand() / (RAND_MAX + 1.0)) - 0.5);
    r_n(2) = r(n,2) + Delta * ((rand() / (RAND_MAX + 1.0)) - 0.5);

    
    double V_tot_r1 = *V;//potenziale in posizione nuova
    
    double dr_mod, dr_mod_n;

    for (int i = 0; i < n; ++i){//ciclo sulla colonna
        for (int k = 0; k < 3; ++k){//ciclo sulle coordinate
            dr_n(i,n,k) = r(i,k) - r_n(k);
            dr_n(i,n,k) -= L * rint(dr_n(i,n,k)/L);//sposto in [-L/2,+L/2]
        }

        dr_mod_n = mod(dr_n,i,n);
        dr_mod = mod(dr,i,n);

        V_tot_r1 += V_LJ(dr_mod_n, L) - V_LJ(dr_mod, L);
    }
    for (int j = n + 1; j < N; ++j){
        for (int k = 0; k < 3; ++k){
            dr_n(n,j,k) = r_n(k) - r(j,k);
            dr_n(n,j,k) -= L * rint(dr_n(n,j,k)/L);//sposto in [-L/2,+L/2]
        }

        dr_mod_n = mod(dr_n,n,j);
        dr_mod = mod(dr,n,j);

        V_tot_r1 += V_LJ(dr_mod_n, L) - V_LJ(dr_mod, L);
    }

    double A = min(1, exp( - (V_tot_r1 - *V) / T));//trovo A
    numero_proposti++;
    
    if(A > (rand() / ( (double)RAND_MAX + 1.0))){
        r.row(n) = r_n;
        numero_accettati++;
        *V = V_tot_r1; 

        *W *= N;
        for (int i = 0; i < n; ++i){//ciclo sulla colonna n
            dr_mod_n = mod(dr_n,i,n);
            dr_mod = mod(dr,i,n);

            dr(i,n,0)=dr_n(i,n,0);//distanze "vecchie per il prossimo ciclo"
            dr(i,n,1)=dr_n(i,n,1);
            dr(i,n,2)=dr_n(i,n,2);

            *W += - dV_LJ(dr_mod_n, L) + dV_LJ(dr_mod, L);
        }
        for (int j = n + 1; j < N; ++j){//ciclo sulla riga n
            dr_mod_n = mod(dr_n,n,j);
            dr_mod = mod(dr,n,j);

            dr(n,j,0)=dr_n(n,j,0);//distanze "vecchie per il prossimo ciclo"
            dr(n,j,1)=dr_n(n,j,1);
            dr(n,j,2)=dr_n(n,j,2);
    
            *W += - dV_LJ(dr_mod_n, L) + dV_LJ(dr_mod, L);
        }
        *W /= N;
    }
}
void calcola_osservabili(rowvec rho, rowvec passi_eq, int caso_min, int N_t, int q, double T_req){
    int caso_max;
    if (caso_min == -1) {
        caso_max = rho.size();
    } else {
        caso_max = caso_min + 1;
    }
    
    //###################################################
    int reticolo   = log2(M);
    gnuplot.open("out/gnuplot.dat");
    gnuplot << caso_min << endl;
    gnuplot.close();

    risultati.open("out/risultati.dat");
    dati.open("out/dati.dat");
    dati_blocking.open("out/dati_blocking.dat");
    int start;
    if(caso_min == -1){
        start = caso_min + 1;
    }
    else{
        start = caso_min;
    }

    mat r(N,3);
    cube dr(N,N,3), dr_n(N,N,3);//li creo all'inizio cosi salvo memoria e tempo

    for (int caso = start; caso < caso_max; ++caso) {

        N_t = passi_eq(caso) + q;//faccio fare 10000 passi dopo l'equilibrazione per avere dei risultati carini

        double L = cbrt(N / rho(caso));
        double Delta;

        Delta = Calcola_Delta(L, rho, caso);

        crea_reticolo(r, L);//creo il reticolo iniziale

        for (int i = 0; i < N; ++i){//creo il cubo di distanze relative alla posizione iniziale
            for (int j = i + 1; j < N; ++j){
                for (int k = 0; k < 3; ++k){
                    dr(i,j,k) = r(i,k) - r(j,k);
                    dr(i,j,k) -= L * rint(dr(i,j,k)/L);//sposto in [-L/2,+L/2]
                }
            }
        }

        double P = 0;
        double E_c = 3./2.*N*T_req;
        double E = 0;
        double V = 0, W = 0;
        double V_m = 0, W_m = 0, P_m = 0, E_m = 0;
        double var_P = 0, var_E = 0;

        // calcolo il primo valore della "pressione" e del potenziale
        for (int i = 0; i < N; ++i){
            for (int j = i + 1; j < N; ++j) {
                double dr_mod = mod(dr,i,j);
                
                V += V_LJ(dr_mod, L);
                W += -dV_LJ(dr_mod, L);
            }
        }
        W /= N;

        for (int i = 0; i < N_t; ++i) {//passi
            if(i > passi_eq(caso) + 1){//calcola le medie
                V_m = V_m * (i - passi_eq(caso) - 1.0);
                W_m = W_m * (i - passi_eq(caso) - 1.0);
                E_m = E_m * (i - passi_eq(caso) - 1.0);

                var_P = var_P * (i - passi_eq(caso) - 1.0);
                var_E = var_E * (i - passi_eq(caso) - 1.0);
            }

            dr_n = dr;
            MRT2(r, &V, &W, N, Delta, T_req, L, dr, dr_n);//metropolis

            E = E_c + V;
            P = (1 + W / (3.0 * T_req));

            if(i > passi_eq(caso)){//tempo di equilibrazione e calcola le medie
                V_m = (V_m + V) / (i - passi_eq(caso));
                W_m = (W_m + W) / (i - passi_eq(caso));
                E_m = (E_m + E) / (i - passi_eq(caso));
                P_m = (1 + W_m / (3.0 * T_req)); 

                var_P = (var_P + (P - P_m) * (P - P_m)) / (i - passi_eq(caso));//calcolo la varianza "ordinaria"
                var_E = (var_E + (E - E_m) * (E - E_m)) / (i - passi_eq(caso));//calcolo la varianza "ordinaria"
            }

            if (caso_min != -1) {
                if(i > passi_eq(caso)){
                    dati << i << "\t" << V_m << "\t" << P_m << "\t" << E_m  << "\t" << V << "\t" << P << "\t" << sqrt(var_P/i) << "\t" << E << "\t" << sqrt(var_E/i) << endl;
                    dati_blocking << P << "\t" << E << endl;
                }
                else{
                    dati << i << "\t" << NAN << "\t" << NAN << "\t" << NAN  << "\t" << V << "\t" << P << "\t" << sqrt(var_P/i) << "\t" << E << "\t" << sqrt(var_E/i) << endl;
                }
            }

            // if(i == (int)(N_t / 2) && caso_min != -1){
            //     cout<<"Sono a met?? dei cicli montecarlo, dai che ce la faccio"<<endl;
            // }
        }

        cout << "Densit?? = "<< rho(caso) << "\t" <<"Pressione = " << P_m << "\t" <<"Energia = " << E_m << endl;
        if (caso_min == -1) {
            dati << rho(caso) << "\t" << P << "\t" << sqrt(var_P) << endl;
        }
        else{
            risultati << rho(caso) << "\t" << P_m << "\t" << sqrt(var_P) << endl;
        }
        cout << "Ho accettato il " << (double)numero_accettati / (double)numero_proposti * 100 << "%. I passi totali erano " << N_t << "\n\n";
    }
    dati_blocking.close();
    dati.close();
    risultati.close();

    if(caso_min != -1){
        stampa_reticolo(r);//stampa il reticolo finale per la gdr
    }
}
// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------
// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------

