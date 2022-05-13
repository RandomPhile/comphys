#include "header.h"
#include "funzioni.h"
#include "es7.h"

/*** variabili globali ***/
//CC, BCC, FCC
int M = 2; //1,2,4
int N = M * pow(5, 3); //numero di particelle
ofstream dati;

struct coppia {
    double rho;
    double sigma;
};

int main() {
    srand(1);//default seed = 1
    double dt      = 0.01;//passo temporale
    double t1      = 10;//durata simulazione

    coppia coppie[] = {//aggiornate con M=2,n=6,dt=0.01,t1=20 (6+6 minuti)
        {.rho = 0.01, .sigma = 1.04249},
        {.rho = 0.1, .sigma = 0.831817},
        {.rho = 0.2, .sigma = 0.687814},
        {.rho = 0.3, .sigma = 0.672185},
        {.rho = 0.4, .sigma = 0.768636},
        {.rho = 0.5, .sigma = 0.921996},
        {.rho = 0.6, .sigma = 1.12905},
        {.rho = 0.7, .sigma = 1.28463},
        {.rho = 0.8, .sigma = 1.41451},
        {.rho = 0.9, .sigma = 1.45276},
        {.rho = 1.0, .sigma = 1.38634},
        {.rho = 1.1, .sigma = 1.39870},
        {.rho = 1.15, .sigma = 1.41261},
        {.rho = 1.2, .sigma = 1.45959}
    };
    int caso = 5;

    //###################################################
    struct vec r[N], v[N], a[N];
    int reticolo   = log2(M);
    double T_req = 1.1;

    double L = cbrt(N / coppie[caso].rho);
    double r_c = L / 2;

    double t = 0;
    int N_t = (t1 - t) / dt;

    crea_reticolo(r, L);
    crea_reticolo(a, 0);
    distr_gauss(v, coppie[caso].sigma);

    vec *dr[N]; vec_2D(dr, N);
    a_LJ(r, a, dr, r_c, L);

    double E = 0, T = 0, P = 0;
    double K = 0, V = 0, W = 0;
    double K_c = 0, V_c = 0, W_c = 0;

    for (int i = 0; i < N_t; ++i) {//tempo
        K = 0, V = 0, W = 0;
        vel_verlet(r, v, a, dt, r_c, L, &K, &V, &W);
        v_cm_0(v);
        t += dt;
    }


    dati.open("dati.dat");
    
    cout << "L = " << L << "\t" << "L/2 = " << L / 2 << "\n" << endl;
    int p = 0;//particella attorno cui calcolo g(r)

    vec rp; double rp_mod[N];
    for (int i = 0; i < N; ++i) {
        rp.x = r[p].x - r[i].x;
        rp.y = r[p].y - r[i].y;
        rp.z = r[p].z - r[i].z;

        rp.x -= L * rint(rp.x / L);
        rp.y -= L * rint(rp.y / L);
        rp.z -= L * rint(rp.z / L);

        rp_mod[i] = rp.mod();
    }
    int N_bins = 15;
    double de_r = 0.5 * L / N_bins;
    cout << "de_r = " << de_r << endl;

    double f;//N(r, dr)
    double g;//g(r)
    for (int j = 0; j < N_bins; ++j) {
        f = 0;
        for (int i = 0; i < N; ++i) {
            //if >= L/2 persi
            if (rp_mod[i] > j * de_r && rp_mod[i] <= (j + 1) * de_r) {
                f++;
            }
        }
        g = f / ((pow((j + 1) * de_r, 3) - pow(j * de_r, 3)) * M_PI * 4.0 / 3.0);
        g /= coppie[caso].rho;
        dati << (2 * j + 1) * de_r / 2 << "\t" << g << "\t" << f << endl;
    }

    dati.close();

    return 0;
}

