#include <iostream>
#include <fstream>
#include <cmath>
#include "eq_diff.h"
#include "funzioni_utili.h"

using namespace std;
ofstream dati;

double f(double t, double x, double y, double *fargs);
double g(double t, double x, double y, double *gargs);
double a(double *r, int indice, double *args);

int M = 1;
int N_mol = M*pow1((double)4, 3);

int main() {
    dati.open("dati.dat");
    double t = 0, h, E, K, T, V;
    double t0 = 0, t1 = 1;

    // int N = 1000;//numero passi temporali
    // h = (t1 - t0) / ( (double) N);//larghezza temporale

    h = 0.001;
    int N = (t1-t0)/h;
    cout<<N<<endl;
    double rho = 1e-2; //fisso la densita del campione da studiare
    double L = cbrt(N_mol / rho);

    double r[3 * N_mol], v[3 * N_mol];//vettori 3N_mol posizioni
    double r_mod, v_mod;

    r_mod = sqrt(r[0] * r[0] + r[0 + 1] * r[0 + 1] + r[0 + 2] * r[0 + 2]);
    v_mod = sqrt(v[0] * v[0] + v[0 + 1] * v[0 + 1] + v[0 + 2] * v[0 + 2]);

    crea_reticolo(N_mol, L, r);//crea il reticolo
    r[0] = r[1] = r[2] = 0;
    stampa_reticolo(N_mol, r);

    gauss_distr(v, 1.1, N_mol);//crea la distribuzione gaussiana di velocita
    set_vcm0(v, N_mol);//impone la velocita media delle particelle a 0
    
    

    double a_prev[3 * N_mol];
    
    compila_matr(r, a_prev, N_mol, fLJ, L);//inizializzo accelerazioni come date dal potenziale al tempo 0
    
    
    E = 0; K = 0; V = 0;
    for (int n = 0; n < N; ++n) {//passi temporali
        K = K * (n + 1); //energia cinetica @t
        V = V * (n + 1); //energia potenziale @t
        
        t = t0 + n * h;//imposto il tempo

        for (int j = 0; j < N_mol; ++j) {//particelle

            if (vel_verlet(t, r, v, h, j, a_prev, fLJ, L)) {printf("ERRORE");}//esegue l'algoritmo di velocity verlet

            //trovo il modulo della particella corrente per trovare il suo contributo alla energia cinetica
            v_mod = sqrt(v[j] * v[j] + v[j + 1] * v[j + 1] + v[j + 2] * v[j + 2]);

            //impongo condizioni di bordo periodiche
            cond_bordo(r, N_mol, L);

            //trovo il valore dell'energia potenziale
            for (int i = 0; i < j; ++i) {//attenzione che è usata la distanza effettiva invece della distanza della prima copia vicina
                double x_ij = r[3 * i] - r[3 * j];
                double y_ij = r[3 * i + 1] - r[3 * j + 1];
                double z_ij = r[3 * i + 2] - r[3 * j + 2];
                V += VLJ(sqrt(x_ij * x_ij + y_ij * y_ij + z_ij * z_ij), L);
            }
            K += 0.5 * v_mod * v_mod;
        }

        K = K / (n + 2);
        V = V / (n + 2);

        E = K + V;
        T = 2.0 * K / (3.0 * N_mol);
        dati << t << "\t" << E << "\t" << K << "\t" << V << "\t" << T << endl;
    }
    dati << "\n\n";

    //stampo il reticolo in posizione finale
    stampa_reticolo(N_mol, r);
    dati.close();
    return 0;
}

double f(double t, double x, double y, double *fargs) {
    return y;
}
double g(double t, double x, double y, double *gargs) {
    return -x;
}
double a(double *r, int indice, double *args) {
    return -r[indice];
}
