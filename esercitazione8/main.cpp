#include <iostream>
#include <cmath>
#include <fstream>

#include "funzioni_utili.h"
#include "distribuzioni.h"

using namespace std;
ofstream dati;

int main() {
    double A=1, B=1;
    double alfa=2, beta=-0.5;
    double parametri_funz[]={alfa, 0, 1};//primo parametro valore moltiplicativo di x, secondo parametro esponente di x davanti alla funzione, terzo parametro valore al quale è elevata l'intera funzione
    double delta=1e-3;
    double mu[4], sigma[4];

    mu[0]=integrale_montecarlo(sin_n, distr_unif, 0, A, delta, parametri_funz);
    parametri_funz[1]=2;
    mu[1]=integrale_montecarlo(sin_n, distr_unif, 0, A, delta, parametri_funz);
    parametri_funz[0]=beta;
    parametri_funz[1]=0;
    mu[2]=2-integrale_montecarlo(exp_n, distr_unif, 0, B, delta, parametri_funz);
    parametri_funz[1]=2;
    mu[3]=16-integrale_montecarlo(exp_n, distr_unif, 0, B, delta, parametri_funz);
    
    for (int i = 0; i < 4; ++i){
        cout<<"il risultato dell'integrale "<<i<<" è "<<mu[i]<<endl;
    }


    //calcolo ora la varianza campionaria
    parametri_funz[2]=2;

    parametri_funz[0]=alfa;
    parametri_funz[1]=0;
    sigma[0]=integrale_montecarlo(sin_n, distr_unif, 0, A, delta, parametri_funz)-pow1(mu[0],2);
    parametri_funz[1]=2;
    parametri_funz[1]=0;
    sigma[1]=integrale_montecarlo(sin_n, distr_unif, 0, A, delta, parametri_funz)-pow1(mu[1],2);
    parametri_funz[0]=2*beta;
    sigma[2]=1-integrale_montecarlo(exp_n, distr_unif, 0, B, delta, parametri_funz)-pow1(mu[2],2);
    parametri_funz[1]=2;
    //sigma[3]=16-integrale_montecarlo_unif(exp_n, 0, B, delta, parametri_funz)-pow1(mu[3],2);

    for (int i = 0; i < 4; ++i){
        cout<<"la varianza campionaria dell'integrale "<<i<<" è "<<sigma[i]<<endl;
    }
}
