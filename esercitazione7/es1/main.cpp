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

<<<<<<< Updated upstream
int N_mol = pow1((double)6, 3);
=======
<<<<<<< HEAD
int N_mol = pow1((double)3, 3);//numero di particelle
=======
int N_mol = pow1((double)6, 3);
>>>>>>> a14900cb9df682e2da9c140b85653d8582dcf285
>>>>>>> Stashed changes

int main() {
    dati.open("dati.dat");
    double t = 0, h, E, K, T, V;
<<<<<<< Updated upstream
    double t0 = 0, t1 = 20;

=======
<<<<<<< HEAD
    double t0 = 0, t1 = 25;

    int N = 1000;//numero di passi temporali
    h = (t1 - t0) / ( (double) N);//larghezza passo temporale
    double rho = 1e-2; //densita del campione da studiare
=======
    double t0 = 0, t1 = 20;

>>>>>>> Stashed changes
    int N = 100;
    h = (t1 - t0) / ( (double) N);
    double rho = 1e-2; //fisso la densita del campione da studiare
>>>>>>> a14900cb9df682e2da9c140b85653d8582dcf285
    double eps, sigma, L, distanza_interaz;
    L = cbrt(N_mol / rho);
    distanza_interaz = L / 2;

<<<<<<< Updated upstream
    double var_ad[] = {eps=1e-6, sigma=1e-1, distanza_interaz, L};
=======
<<<<<<< HEAD
    double var_ad[] = {eps = 1e-6, sigma = 1e-6, distanza_interaz, L};
=======
    double var_ad[] = {eps=1e-6, sigma=1e-1, distanza_interaz, L};
>>>>>>> a14900cb9df682e2da9c140b85653d8582dcf285
>>>>>>> Stashed changes

    double r[3 * N_mol], v[3 * N_mol];//vettori 3N dimensionali posizione e velocità
    double r_mod, v_mod;//scalari modulo

    r_mod = sqrt(r[0] * r[0] + r[0 + 1] * r[0 + 1] + r[0 + 2] * r[0 + 2]);
    v_mod = sqrt(v[0] * v[0] + v[0 + 1] * v[0 + 1] + v[0 + 2] * v[0 + 2]);

<<<<<<< Updated upstream
    crea_reticolo(N_mol, L, r);//crea il reticolo
    
    gauss_distr(v, 1.1, N_mol);//crea la distribuzione gaussiana di velocita
    set_vcm0(v, N_mol);//impone la velocita media delle particelle a 0

    double a_prev[3 * N_mol];
=======
<<<<<<< HEAD
    crea_reticolo(N_mol, L, r);//dispone i valori di posizione lungo un reticolo
    stampa_reticolo(N_mol, r);//salva nel file dati i valori delle posizioni

    gauss_distr(v, 1.2, N_mol);//aggiorno le velocità secondo una distribuzione gaussiana
    set_vcm0(v, N_mol);

    double a_prev[3 * N_mol];
    compila_matr(r, a_prev, N_mol, fLJ, var_ad);//inizializzo accelerazioni come date dal potenziale

=======
    crea_reticolo(N_mol, L, r);//crea il reticolo
    
    gauss_distr(v, 1.1, N_mol);//crea la distribuzione gaussiana di velocita
    set_vcm0(v, N_mol);//impone la velocita media delle particelle a 0

    double a_prev[3 * N_mol];
>>>>>>> Stashed changes
    compila_matr(r, a_prev, N_mol, fLJ, var_ad);//inizializzo accelerazioni come date dal potenziale al tempo 0
    
>>>>>>> a14900cb9df682e2da9c140b85653d8582dcf285
    E = 0; K = 0; V = 0;
    for (int n = 0; n < N; ++n) {//passi temporali
        K = K * (n + 1); //energia cinetica @t
        V = V * (n + 1); //energia potenziale @t
        T = 0;//temp @t
<<<<<<< Updated upstream
        
        t = t0 + n * h;//imposto il tempo
        
        for (int j = 0; j < N_mol; ++j) {//particelle
            
            if (vel_verlet(t, r, v, h, j, a_prev, fLJ, var_ad)) {printf("ERRORE");}//esegue l'algoritmo di velocity verlet
            
            //trovo il modulo della particella corrente per trovare il suo contributo alla energia cinetica
            v_mod = sqrt(v[j] * v[j] + v[j + 1] * v[j + 1] + v[j + 2] * v[j + 2]);
=======
<<<<<<< HEAD
        //           r_x           r_y                   r_z
        r_mod = sqrt(r[0] * r[0] + r[0 + 1] * r[0 + 1] + r[0 + 2] * r[0 + 2]);
        v_mod = sqrt(v[0] * v[0] + v[0 + 1] * v[0 + 1] + v[0 + 2] * v[0 + 2]);

        for (int j = 0; j < N_mol; ++j) {//per N particelle j
            t = t0 + n * h;//incremento tempo

            //aggiorno posizione e velocità con vel_verlet:
            if (vel_verlet(t, r, v, h, j, a_prev, fLJ, var_ad)) {printf("ERRORE");}

            v_mod = sqrt(v[j] * v[j] + v[j + 1] * v[j + 1] + v[j + 2] * v[j + 2]);


            cond_bordo(r, N_mol, L);//aggiorno la posizione spostando le particelle dentro alla scatola
            for (int i = 0; i < j; ++i) {//attenzione che è usata la distanza effettiva invece della distanza della prima copia vicina
                double x_ij = r[3 * i] - r[3 * j];
                double y_ij = r[3 * i + 1] - r[3 * j + 1];
                double z_ij = r[3 * i + 2] - r[3 * j + 2];
                V += VLJ(sqrt(x_ij * x_ij + y_ij * y_ij + z_ij * z_ij), var_ad);
            }

=======
        
        t = t0 + n * h;//imposto il tempo
        
        for (int j = 0; j < N_mol; ++j) {//particelle
            
            if (vel_verlet(t, r, v, h, j, a_prev, fLJ, var_ad)) {printf("ERRORE");}//esegue l'algoritmo di velocity verlet
            
            //trovo il modulo della particella corrente per trovare il suo contributo alla energia cinetica
            v_mod = sqrt(v[j] * v[j] + v[j + 1] * v[j + 1] + v[j + 2] * v[j + 2]);
>>>>>>> Stashed changes
            
            //impongo condizioni di bordo periodiche
            cond_bordo(r, N_mol, L);
            
            //trovo il valore dell'energia potenziale
            for (int i = 0; i < j; ++i) {//attenzione che è usata la distanza effettiva invece della distanza della prima copia vicina 
                double x_ij=r[3 * i] - r[3 * j];
                double y_ij=r[3 * i + 1] - r[3 * j + 1];
                double z_ij=r[3 * i + 2] - r[3 * j + 2];
                V += VLJ(sqrt(x_ij * x_ij + y_ij * y_ij + z_ij * z_ij), var_ad);
            }
<<<<<<< Updated upstream
=======
>>>>>>> a14900cb9df682e2da9c140b85653d8582dcf285
>>>>>>> Stashed changes
            K += 0.5 * v_mod * v_mod;
            
        }
        
        K = K / (n + 2);
        V = V / (n + 2);

        E = K + V;
        T = 2.0 * K / (3.0 * N_mol);
        dati << t  << "\t" << E << "\t" << K << "\t" << V << "\t" << T << endl;
        //cout << t  << "\t" << E << "\t" << K << "\t" << V << "\t" << T << endl;
    }
<<<<<<< HEAD
    dati << "\n\n";
    stampa_reticolo(N_mol, r);
=======
    dati<<"\n\n";
        
    //stampo il reticolo in posizione finale
    stampa_reticolo(N_mol, L, r);
>>>>>>> a14900cb9df682e2da9c140b85653d8582dcf285
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

