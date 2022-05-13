#include "header.h"
#include "funzioni.h"
#include "es6.h"

/*** variabili globali ***/
int N = 100; //numero di particelle (meglio >100)
ofstream dati;

int main() {
    srand(4);//default seed = 1
    struct vec r[N], v[N], a[N];

    double t = 0, t1 = 25;
    double dt = 0.01;
    int N_t = (t1 - t) / dt;

    for (int i = 0; i < N; ++i) {
        r[i].uguale(0);
        v[i].uguale(0);
        a[i].uguale(0);
    }
    double sigma = 1;
    distr_gauss(v, sigma);
    //distr_gauss(r, sigma);//parte 6.2b
    /*
    scegliendo come posizione iniziale l'origine
    usando sigma = 1
    T si equilibra a ~0.5 (un po' meno)
    
    distribuendo invece la pos iniziale come le velocità
    bisogna usare sigma = 0.71
    per ottenere la stessa T_eq
    */
    
    v_cm_0(v);


    dati.open("dati.dat");
    double K = 0, V = 0, E = 0, T = 0;
    double K_c = 0;

    for (int i = 0; i < N_t; ++i) {//tempo
        K_c = K_c * (i + 1.0);

        vel_verlet(r, v, a, dt, &K, &V);

        K_c = (K_c + K) / (i + 2.0);
        //stampa_stato(r,v,a);
        T = 2.0 * K_c / (3.0 * N);
        E = K + V;
        dati << t << "\t" << K << "\t" << V << "\t" << E << "\t" << T << endl;
        t += dt;
    }

    dati.close();
    return 0;
}


