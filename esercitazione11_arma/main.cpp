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
int M = 4; //1,2,4
int N = M * pow(8, 3); //numero di particelle

int numero_proposti=0;
int numero_accettati=0;
ofstream dati, gnuplot, coord;

int main() {
    srand(1);//default seed = 1
    int N_t;//numero passi simulazione, vedi rigo 61 per il valore in base al tempo di equilibrazione stimato
    int N_b=20;//numero bin
    double T_req = 1.1;//temperatura adimensionale

    rowvec rho={0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
    rowvec passi_eq={400, 550, 1200, 1000, 2200, 2500, 2500, 2800, 3700, 3700, 3700, 3700, 3700, 3700, 3600, 2200};//tempi di equilibrazione con delta preso per avere circa 50%

    int caso_min = 0;//mettere -1 per avere P(rho)
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

    coord.open("coord.dat");
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

        N_t=passi_eq(caso)+1e4;//faccio fare 7000 passi dopo l'equilibrazione per avere dei risultati carini
        //N_t=3e4;

        double L = cbrt(N / rho(caso));
        double r_c = L / 2;
        double Delta;

        if(caso>=8){
            Delta=L/(60*rho(caso)*rho(caso));//scelgo un delta che mi dia circa 50% di accettazione
            if(caso!=-1){
                cout<<"Delta scalato sulla scatola: "<<Delta/L<<endl;
            }
        }
        else if(caso<8 && caso>4){
            Delta=L/(50*rho(caso));//scelgo un delta che mi dia circa 50% di accettazione
            if(caso!=-1){
                cout<<"Delta scalato sulla scatola: "<<Delta/L<<endl;
            }
        }
        else{
            Delta=L/(70*rho(caso));//scelgo un delta che mi dia circa 50% di accettazione
            if(caso!=-1){
                cout<<"Delta scalato sulla scatola: "<<Delta/L<<endl;
            }
        }
        
        

        double P = 0, var_P = 0;
        double V = 0, W = 0;
        double V_m = 0, W_m = 0, P_m = 0;

        calcola_coord_oss();

        cout << "Densità = "<< rho(caso) << "\t" <<"Pressione = " << P << endl;
        if (caso_min == -1) {
            dati << rho(caso) << "\t" << P << endl;
        }
        else{//calcolo della gdr
            gdr_funz(r, L, rho(caso), N_b);
        }
        

        cout <<"Ho accettato il " <<(double)numero_accettati/(double)numero_proposti*100 << "%. I passi totali erano "<<N_t<<"\n\n";
    }
    
    dati.close();
    coord.close();


    return 0;
}

