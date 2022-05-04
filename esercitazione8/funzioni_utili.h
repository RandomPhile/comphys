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
double integrale_montecarlo_unif(double (*f)(double,double*), double inizio, double fine, double errore, double *args){
    int N=rint(pow1(1/errore,2)+0.5);
    //cout<<"N= "<<N<<endl;
    double ris=0;
    for (int i = 1; i <= N; ++i){
        double x=(rand()/((double)RAND_MAX+1.0))*abs(fine-inizio)+inizio;//casuale tra zero e uno uniforme
        ris+=f(x,args);
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

#endif /* funzioni_utili_h */