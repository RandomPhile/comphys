//
//  funzioni_utili.h
//  es1
//  Created by Mattia Lupi on 18/04/22.
//

#ifndef funzioni_utili_h
#define funzioni_utili_h

using namespace std;

extern ofstream dati;

double pow1(double base, int esp) {//funzione esponenziale creata per non usare pow
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}
void plot_hist(double *v, int N_mol) {
    //usare N grande
    double v_max = 0;
    double v_min = 1e50;
    int N_v = 1 + 3.322 * log(N_mol); //Sturge's Rule
    for (int i = 0; i < N_mol; ++i){
        if(v[i]<v_min){
            v_min=v[i];
        }
        if(v[i]>v_max){
            v_max=v[i];
        }
    }
    double v_step = (v_max - v_min) / N_v;
    double bins;
    double f_v;
    for (int j = 0; j < N_v; ++j) {
        bins = v_min + v_step * j;
        f_v = 0;
        for (int i = 0; i < N_mol; ++i) {
            if (v[i] >= bins && v[i] < (bins + v_step)) {
                f_v++;
            }
        }
        dati << bins << "\t" << f_v << endl;
    }
    dati << "\n\n";
}

double integrale_montecarlo(double (*f)(double, double*), void (*distr)(double*, int, double*), double errore, double *args_f, double *args_d){//versione distrib uniforme
    //args_d: 0 inizio, 1 fine, 2 sigma, 3 mu
    
    int N=rint(pow1(1/errore,2)+0.5);//
    double ris=0;
    double x[N];
    distr(x, N, args_d);
    for (int i = 1; i <= N; ++i){
        ris+=f(x[i],args_f);
    }
    ris=ris*(args_d[1]-args_d[0])/N;
    return ris;
}
double integrale_montecarlo(double (*f)(double, double*), double *args_f,
                            void (*distr)(double*,int, double*), double *args_d,
                            double errore,
                            double (*g)(double, double*), double *args_g){//versione g(x) generica
    //args_d: 0 inizio, 1 fine, 2 sigma, 3 mu
    
    int N=rint(  pow1(args_d[2]/errore,2)   +0.5   );//
    double ris=0;
    double x[N];
    distr(x, N, args_d);
    //plot_hist(x,N);
    for (int i = 0; i <  N; ++i){
        ris+=f(x[i],args_f)/g(x[i],args_g);
    }
    ris=ris/N;
    return ris;
}


double sin_n(double x, double *par){//0 alfa interno al seno, 1 l'esponente di x, 2 il valore di cui è elevato il tutto
    double ris=pow1(sin(x*par[0]),2)*pow1(x, par[1]);
    return pow1(ris, par[2]);
}
double exp_n(double x, double *par){//0 alfa interno all'esponenziale, 1 l'esponente di x, 2 il valore di cui è elevato il tutto
    double ris=exp(x*par[0])*pow1(x, par[1]);
    return pow1(ris, par[2]);
}
double gauss(double x, double *par){//primo parametro sigma, secondo parametro mu
    double ris=1/(sqrt(2*M_PI)*par[0])*exp(-(x-par[1])*(x-par[1])/(2*par[0]*par[0]));
    return ris;
}
double g_prof(double x, double *par){
    if(x>=3){
        return 1;
    }
    else{
        return 0;
    }
}

#endif /* funzioni_utili_h */
