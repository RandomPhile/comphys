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
int M = 1; //1,2,4
int N = M * pow(7, 3); //numero di particelle

int numero_proposti=0;
int numero_accettati=0;
ofstream dati, gnuplot, risultati;











// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------
// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------

double mod(cube &r, double riga, double colonna);
double V_LJ(double r, double L);
void MRT2(mat &r, double *V, double *W, int N, double Delta, double T, double L, cube &dr);

// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------
// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------
















// -------------------------------------------MAIN---------------------------------MAIN-----------------------------------MAIN-----------------------------------------------
// -------------------------------------------MAIN---------------------------------MAIN-----------------------------------MAIN-----------------------------------------------

int main() {
    srand(1);//default seed = 1
    int N_t;//numero passi simulazione, vedi rigo 61 per il valore in base al tempo di equilibrazione stimato
    int N_b=20;//numero bin
    double T_req = 1.1;//temperatura adimensionale

    rowvec rho={0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
    rowvec passi_eq={400, 550, 1200, 1000, 2200, 2500, 2500, 2800, 3700, 3700, 3700, 3700, 3700, 3700, 3600, 2200};//tempi di equilibrazione con delta preso per avere circa 50%



    int caso_min = 8;//mettere -1 per avere P(rho)



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
    gnuplot.open("gnuplot.dat");
    gnuplot << caso_min << endl;
    gnuplot.close();

    int q=1e5;//numero punti in piu rispetto al tempo di equilibrazione


    // risultati.open("risultati.dat");
    // dati.open("dati.dat");
    // int start;
    // if(caso_min==-1){
    //     start=caso_min+1;
    // }
    // else{
    //     start=caso_min;
    // }
    // for (int caso = start; caso < caso_max; ++caso) {

    //     N_t=passi_eq(caso)+q;//faccio fare 10000 passi dopo l'equilibrazione per avere dei risultati carini

    //     double L = cbrt(N / rho(caso));
    //     double Delta;

    //     Delta=Calcola_Delta(L, rho, caso);

    //     crea_reticolo(r, L);//creo il reticolo iniziale

    //     for (int i = 0; i < N; ++i){//creo il cubo di distanze relative alla posizione iniziale
    //         for (int j = i + 1; j < N; ++j){
    //             for (int k = 0; k < 3; ++k){
    //                 dr(i,j,k) = r(i,k) - r(j,k);
    //                 dr(i,j,k) -= L * rint(dr(i,j,k)/L);//sposto in [-L/2,+L/2]
    //             }
    //         }
    //     }

    //     double P = 0, sigma_P;
    //     double V = 0, W = 0;
    //     double V_m = 0, W_m = 0, P_m = 0;
    //     double var_P = 0;

    //     for (int i = 0; i < N_t; ++i) {//passi
    //         if(i>passi_eq(caso)){//calcola le medie
    //             V_m = V_m * (i - passi_eq(caso));
    //             W_m = W_m * (i - passi_eq(caso));
    //             var_P = var_P * (i - passi_eq(caso));
    //         }

    //         MRT2(r, &V, &W, N, Delta, T_req, L, dr);//metropolis
    //         // cout<<V<<endl;

    //         P = rho(caso) * (1 + W / (3.0 * T_req));
    //         if(i>passi_eq(caso)){//tempo di equilibrazione e calcola le medie
    //             V_m = (V_m + V) / (i - passi_eq(caso) + 1.0);
    //             W_m = (W_m + W) / (i - passi_eq(caso) + 1.0);

                
    //             P_m = rho(caso) * (1 + W_m / (3.0 * T_req)); //P su rho*k_B*T_req
    //             var_P = (var_P + (P - P_m) * (P - P_m)) / (i - passi_eq(caso) + 1.0);//calcolo la varianza "ordinaria"
    //         }

    //         if (caso_min != -1) {
    //             dati << i << "\t" << V_m << "\t" << P_m  << "\t" << V << "\t" << P << "\t" << sqrt(var_P) << endl;
    //         }

    //         if(i==(int)(N_t/2) && caso_min!=-1){
    //             cout<<"Sono a metà dei cicli montecarlo, dai che ce la faccio"<<endl;
    //         }
    //     }

    //     cout << "Densità = "<< rho(caso) << "\t" <<"Pressione = " << P_m << endl;
    //     if (caso_min == -1) {
    //         dati << rho(caso) << "\t" << P << "\t" << sqrt(var_P) << endl;
    //     }
    //     else{//calcolo della gdr
    //         risultati << rho(caso) << "\t" << P_m << "\t" << sqrt(var_P) << endl;
    //        // gdr_funz(r, L, rho(caso), N_b);
    //     }
    //     cout <<"Ho accettato il " <<(double)numero_accettati/(double)numero_proposti*100 << "%. I passi totali erano "<<N_t<<"\n\n";
    // }
    // // plot_pressioni(caso_min);
    // dati.close();
    // risultati.close();

    blocking(passi_eq(caso_min)+q);

    // blocking_plot();
    return 0;
}

// -------------------------------------------MAIN---------------------------------MAIN-----------------------------------MAIN-----------------------------------------------
// -------------------------------------------MAIN---------------------------------MAIN-----------------------------------MAIN-----------------------------------------------


























// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------
// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------
double mod(cube &r, double riga, double colonna){//calcolo il modulo della posizione relativa delle particelle i e j
    double mod= sqrt(pow1(r(riga,colonna,0),2)+pow1(r(riga,colonna,1),2)+pow1(r(riga,colonna,2),2));
    return mod;
}
double V_LJ(double r, double L) {
    //potenziale va a zero in modo continuo sul bordo della scatola
    if (r < L / 2 && r != 0) {
        double VL_2 = 4 * (pow1(2 / L, 12) - pow1(2 / L, 6));
        //double VpL_2 = 24 * (pow1(1 / r, 8) - 2 * pow1(1 / r, 14));
        return 4 * (pow1(1 / r, 12) - pow1(1 / r, 6)) - VL_2;
    } else {
        return 0;
    }
}
void MRT2(mat &r, double *V, double *W, int N, double Delta, double T, double L, cube &dr){//N=N_mol
    cube dr_n(N, N, 3);
    rowvec r_n(3); 

    double s= rand() / (RAND_MAX + 1.0);//eseguo prima il rand perche da problemi senno
    int n = (int)rint((N-1) *s);//trovo la molecola che viene modificata da MTR2

    r_n(0) = r(n,0) + Delta * (rand() / (RAND_MAX + 1.0) - 0.5);//modifico le posizioni della particella n
    r_n(1) = r(n,1) + Delta * (rand() / (RAND_MAX + 1.0) - 0.5);
    r_n(2) = r(n,2) + Delta * (rand() / (RAND_MAX + 1.0) - 0.5);
    
    double V_tot_r1= *V;//potenziale in posizione nuova
    double V_tot_r0= *V;//potenziale in posizione vecchia
    double dr_mod, dr_mod_n;
    
    dr_n=dr;

    for (int i = 0; i < n; ++i){//ciclo sulla colonna
        for (int k = 0; k < 3; ++k){//ciclo sulle coordinate
            dr_n(i,n,k) = r(i,k) - r_n(k);
            dr_n(i,n,k) -= L * rint(dr_n(i,n,k)/L);//sposto in [-L/2,+L/2]
        }

        dr_mod_n = mod(dr_n,i,n);
        dr_mod = mod(dr,i,n);

        V_tot_r1+=V_LJ(dr_mod_n, L) - V_LJ(dr_mod, L);
    }
    for (int j = n + 1; j < N; ++j){
        for (int k = 0; k < 3; ++k){
            dr_n(n,j,k) = r_n(k) - r(j,k);
            dr_n(n,j,k) -= L * rint(dr_n(n,j,k)/L);//sposto in [-L/2,+L/2]
        }

        dr_mod_n = mod(dr_n,n,j);
        dr_mod = mod(dr,n,j);

        V_tot_r1+=V_LJ(dr_mod_n, L) - V_LJ(dr_mod, L);
    }

    double A=min(1,exp(-(V_tot_r1-V_tot_r0)/T));//trovo A
    numero_proposti++;
    
    if(A>(rand()/((double)RAND_MAX+1.0))){
        r.row(n)=r_n;
        numero_accettati++;
        *V = V_tot_r1; *W = 0;
        for (int i = 0; i < N; ++i){
            for (int j = i + 1; j < N; ++j) {
                dr_mod=mod(dr_n,i,j);
                dr(i,j,0)=dr_n(i,j,0);
                dr(i,j,1)=dr_n(i,j,1);
                dr(i,j,2)=dr_n(i,j,2);
                
                if (dr_mod < L/2) {    
                    //dV_dr * r
                    *W -= 24 * (pow1(1 / dr_mod, 6) - 2 * pow1(1 / dr_mod, 12)); 
                }
            }
        }
        *W/=N;
    }
}
// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------
// -----------------------------------------FUNZIONI----------------------------FUNZIONI------------------------------FUNZIONI-----------------------------------------------

