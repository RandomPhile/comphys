//
//  funzioni_utili.h
//  es1
//  Created by Mattia Lupi on 18/04/22.
//

#ifndef funzioni_utili_h
#define funzioni_utili_h
using namespace std;

double pow1(double base, int esp) {//funzione esponenziale creata per non usare pow
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}
double integrale_montecarlo(double (*f)(double, double*), void (*distr)(double*, int, double, double), double inizio, double fine, double errore, double *args_f){
    int N=rint(pow1(1/errore,2)+0.5);//double *x, int L, double a, double b
    //cout<<"N= "<<N<<endl;
    double ris=0;
    double x[N];
    distr(x, N, inizio, fine);
    for (int i = 1; i <= N; ++i){
        ris+=f(x[i],args_f);
    }
    ris=ris*(fine-inizio)/N;
    return ris;
}

double sin_n(double x, double *par){//primo parametro il coefficiente alfa interno al seno, secondo parametro l'esponente di x, terzo parametro il valore di cui è elevato il tutto
    double ris=pow1(sin(x*par[0]),2)*pow1(x, par[1]);
    return pow1(ris, par[2]);
}
double exp_n(double x, double *par){//primo parametro il coefficiente alfa interno all'esponenziale, secondo parametro l'esponente di x, terzo parametro il valore di cui è elevato il tutto
    double ris=exp(x*par[0])*pow1(x, par[1]);
    return pow1(ris, par[2]);
}
double gauss(double x, double *par){//primo parametro sigma, secondo parametro mu
    double ris=exp((x-par[1])*(x-par[1])/(2*par[0]*par[0]));
    return ris;
}

#endif /* funzioni_utili_h */
