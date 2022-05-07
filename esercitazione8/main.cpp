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
    double parametri_funz[]={alfa, 0, 1};//0 valore moltiplicativo di x, 1 esponente di x davanti alla funzione, 2 valore al quale è elevata l'intera funzione
    double par_distr[]={0, A, 1, A/2};// 0 inizio, 1 fine, 2 sigma, 3 mu
    double delta=1e-3;
    double mu[4], sigma[4];

    mu[0]=integrale_montecarlo(sin_n, distr_unif, delta, parametri_funz, par_distr);
    parametri_funz[1]=2;
    mu[1]=integrale_montecarlo(sin_n, distr_unif, delta, parametri_funz, par_distr);
    parametri_funz[0]=beta;
    parametri_funz[1]=0;
    mu[2]=2-integrale_montecarlo(exp_n, distr_unif, delta, parametri_funz, par_distr);
    parametri_funz[1]=2;
    mu[3]=16-integrale_montecarlo(exp_n, distr_unif, delta, parametri_funz, par_distr);
    
    for (int i = 0; i < 4; ++i){
        cout<<"il risultato dell'integrale "<<i<<" è "<<mu[i]<<endl;
    }


    //calcolo ora la varianza campionaria
    parametri_funz[2]=2;

    parametri_funz[0]=alfa;
    parametri_funz[1]=0;
    sigma[0]=integrale_montecarlo(sin_n, distr_unif, delta, parametri_funz, par_distr)-pow1(mu[0],2);
    parametri_funz[1]=2;
    parametri_funz[1]=0;
    sigma[1]=integrale_montecarlo(sin_n, distr_unif, delta, parametri_funz, par_distr)-pow1(mu[1],2);
    parametri_funz[0]=2*beta;
    sigma[2]=1-integrale_montecarlo(exp_n, distr_unif, delta, parametri_funz, par_distr)-pow1(mu[2],2);
    parametri_funz[1]=2;
    //sigma[3]=16-integrale_montecarlo_unif(exp_n, delta, parametri_funz, par_distr)-pow1(mu[3],2);

    for (int i = 0; i < 4; ++i){
        cout<<"la varianza campionaria dell'integrale "<<i<<" è "<<sigma[i]<<endl;
    }
    
    par_distr[0]=3;//inizio
    par_distr[1]=5;//fine
    par_distr[2]=0.5;//sigma distr
    par_distr[3]=3.7;//mu distr
    parametri_funz[0]=1; //sigma funz
    parametri_funz[1]=0; //mu funz
    double par_g[]={0.5, 3.7};//sigma e mu distr
    
    cout<<"\n\n"<<"Integrale con gaussiana "<<integrale_montecarlo(gauss, parametri_funz, distr_gauss, par_distr, delta, gauss, par_g)<<endl;
    cout<<"Integrale con gaussiana + specchiamento "<<integrale_montecarlo(gauss, parametri_funz, distr_gauss_mezza, par_distr, delta, gauss, par_g)<<endl;
    
}
