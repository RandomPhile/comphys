#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

// #include "funzioni_utili.h"
//#include "distribuzioni.h"

using namespace std;
using namespace arma;

ofstream dati;

void plot_hist(mat v, int N_mol);
double pow1(double base, int esp);
void distr_gauss(mat x, int L, double *args);

int main() {
    dati.open("dati.dat");
    int N=pow1(6,3);
    int M=1;
    double args[]={0,1,0,1};
    N*=M;
    Mat<double> v(3*N,1);
    v.randn();
    cout<<"v=\n"<<v.randn()<<endl;
    plot_hist(v,3*N);
    // distr_gauss(v, 3*N, args);
    // plot_hist(v,3*N);
    dati.close();
    return 0;
}
void plot_hist(mat v, int N_mol) {
    //usare N grande
    double v_max = 0;
    double v_min = 1e50;
    int N_v = 1 + 3.322 * log(N_mol); //Sturge's Rule
    for (int i = 0; i < N_mol; ++i){
        if(v(i,1)<v_min){
            v_min=v(i,1);
        }
        if(v(i,1)>v_max){
            v_max=v(i,1);
        }
    }
    double v_step = (v_max - v_min) / N_v;
    double bins;
    double f_v;
    for (int j = 0; j < N_v; ++j) {
        bins = v_min + v_step * j;
        f_v = 0;
        for (int i = 0; i < N_mol; ++i) {
            if (v(i,1) >= bins && v(i,1) < (bins + v_step)) {
                f_v++;
            }
        }
        dati << bins << "\t" << f_v << endl;
    }
    dati << "\n\n";
}
double pow1(double base, int esp) {//funzione esponenziale creata per non usare pow
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}
void distr_gauss(mat x, int L, double *args) {
    //Formule di Box-Muller
    double x1, x2;
    for (int i = 0; i < L; i += 2) {
        x1 = rand() / (RAND_MAX + 1.0);
        x2 = rand() / (RAND_MAX + 1.0);

        x(i,1)     = args[2] * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1) + args[3];
        x(i + 1,1) = args[2] * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1) + args[3];
    }
}