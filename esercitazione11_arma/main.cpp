#include <iostream>
#include <fstream>
#include <cmath>

#define ARMA_NO_DEBUG

#include <armadillo>
#include <iomanip>
#include <algorithm>

#include "es11.h"

/* NOTE:
-se K o V esplodono molto probabilmente dt è troppo grande
-per stimare sigma t.c. T=1.1, uso N grande ()
-M=1, caso=2:  non esiste un sigma per far convergere T a T_req=1.1
*/

/*** variabili globali ***/
//CC, BCC, FCC
int M = 2; //1,2,4
int N = M * pow(6, 3); //numero di particelle

int numero_proposti=0;
int numero_accettati=0;
ofstream dati, gnuplot, risultati;

int main() {
    srand(1);//default seed = 1
    int N_t = 1e5;//numero passi simulazione
    int N_b=30;//numero bin
    double T_req = 1.1;//temperatura adimensionale

    rowvec rho={0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4};

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
        double L = cbrt(N / rho(caso));
        double r_c = L / 2;

        double Delta=L/(40*rho(caso));//scelgo un delta che mi dia circa 50% di accettazione
        
        crea_reticolo(r, L);//creo il reticolo iniziale
        
        for (int i = 0; i < N; ++i){//creo il cubo di distanze relative alla posizione iniziale
            for (int j = i + 1; j < N; ++j){
                for (int k = 0; k < 3; ++k){
                    dr(i,j,k) = r(i,k) - r(j,k);
                    dr(i,j,k) -= L * rint(dr(i,j,k)/L);//sposto in [-L/2,+L/2]
                }
            }
        }

        double E = 0, T = 0, P = 0;
        double K = 0, V = 0, W = 0;
        double K_c = 0, V_m = 0, W_m = 0, P_m = 0;

        for (int i = 0; i < N_t; ++i) {//tempo
            V_m = V_m * (i + 1.0);
            W_m = W_m * (i + 1.0);
            
            MRT2(r, &V, &W, N, Delta, T_req, L, r_c, dr);
            
            V_m = (V_m + V) / (i + 2.0);
            W_m = (W_m + W) / (i + 2.0);

            P = (1 + W / (3.0 * T_req));
            P_m = (1 + W_m / (3.0 * T_req)); //P su rho*k_B*T_req
            // P = rho(caso) * (1 + W_c / (3.0 * T_req)); //P su k_B*T_req

            if (caso_min != -1) {
                dati << i << "\t" << V_m << "\t" << P_m  << "\t" << V << "\t" << P << endl;
            }
            if(i==(int)(N_t/2) && caso_min!=-1){
                cout<<"Sono a metà dei cicli montecarlo, dai che ce la faccio"<<endl;
            }
        }
        
        if (caso_min == -1) {
            dati << rho(caso) << "\t" << P << endl;
            cout << "Densità = "<< rho(caso) << "\t" <<"Pressione = " << P << endl;
        }
        else{//calcolo della gdr
            cout<<"Ora calcolo la g(r), abbi ancora un po' di pazienza"<<endl;
            mat gdr(N_b, 2);//in una la gdr e nell'altra colonna la distanza dalla part centrale
            
            dati<<"\n\n";
            double L_cella= L / cbrt(N / M);
            gdr_funz(r, L, rho(caso), gdr, N_b);
            for (int i = 0; i < N_b; ++i){
                dati << gdr(i,1)<< "\t" << gdr(i,0) << endl;
            }
            dati<<"\n\n";
            for (int i = 0; i < N; ++i){
                dati<<r(i,0)/L<<"\t"<<r(i,1)/L<<"\t"<<r(i,2)/L<<endl;//normalizzo sulla scatola per presentare il risultato
            }
        }
        

        cout <<"Ho accettato il " <<(double)numero_accettati/(double)numero_proposti*100 << "%\n\n";
    }
    
    dati.close();
    risultati.close();


    return 0;
}

