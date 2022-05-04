#include <iostream>
#include <cmath>
#include <fstream>

#include "funzioni_utili.h"

using namespace std;
ofstream dati;

int main() {
    double A=1, B=1;
    double alfa=2, beta=-0.5;
    double parametri[]={alfa, 0, 1};
    double delta=1e-3;
    double mu[4], sigma[4];

    mu[0]=integrale_montecarlo_unif(sin_n, 0, A, delta, parametri);
    parametri[1]=2;
    mu[1]=integrale_montecarlo_unif(sin_n, 0, A, delta, parametri);
    parametri[0]=beta;
    mu[2]=2-integrale_montecarlo_unif(exp_n, 0, B, delta, parametri);
    parametri[1]=2;
    mu[3]=16-integrale_montecarlo_unif(exp_n, 0, B, delta, parametri);
    
    for (int i = 0; i < 4; ++i){
        cout<<"il risultato dell'integrale "<<i<<" è "<<mu[i]<<endl;
    }


    //calcolo ora la varianza campionaria
    parametri[2]=2;

    parametri[0]=alfa;
    parametri[1]=0;
    sigma[0]=integrale_montecarlo_unif(sin_n, 0, A, delta, parametri)-pow1(mu[0],2);
    parametri[1]=2;
    sigma[1]=integrale_montecarlo_unif(sin_n, 0, A, delta, parametri)-pow1(mu[1],2);
    //parametri[0]=beta;
    //sigma[2]=2-integrale_montecarlo_unif(exp_n, 0, B, delta, parametri)-pow1(mu[2],2);
    //parametri[1]=2;
    //sigma[3]=16-integrale_montecarlo_unif(exp_n, 0, B, delta, parametri)-pow1(mu[3],2);

    for (int i = 0; i < 4; ++i){
        cout<<"la varianza campionaria dell'integrale "<<i<<" è "<<sigma[i]<<endl;
    }
}
