#ifndef funzioni_h
#define funzioni_h


#include "header.h"
struct vec {
    double x; double y; double z;

    double mod() {
        return sqrt(x * x + y * y + z * z);
    }
    void uguale(double val) {
        x = val; y = val; z = val;
    }
    void piu(double val) {
        x += val; y += val; z += val;
    }
    void per(double val) {
        x *= val; y *= val; z *= val;
    }
};
void vec_2D (vec *dr[], int dim) {
    for (int i = 0; i < dim; ++i)
        dr[i] = new vec[dim];
}
void stampa_stato(vec *r, vec *v, vec *a) {
    for (int i = 0; i < N; ++i) {
        cout << "i = " << i << endl;
        cout << "rx = " << r[i].x << endl;
        cout << "vx = " << v[i].x << endl;
        cout << "ax = " << a[i].x << "\n" << endl;

        cout << "ry = " << r[i].y << endl;
        cout << "vy = " << v[i].y << endl;
        cout << "ay = " << a[i].y << "\n" << endl;

        cout << "rz = " << r[i].z << endl;
        cout << "vz = " << v[i].z << endl;
        cout << "az = " << a[i].z << "\n" << endl;
        cout << "\n" << endl;
    }
}
void stampa_costanti(int M, int N, int caso, double *rho, double L, double t1, double dt, int N_t, int N_step, int skip, double pausa, double sigma) {
    cout << "M      = " << M << endl;
    cout << "N      = " << N << endl;
    cout << "caso   = " << caso << endl;
    cout << "rho    = " << rho[caso] << endl;
    cout << "L      = " << L << endl;
    cout << "t1     = " << t1 << endl;
    cout << "dt     = " << dt << endl;
    cout << "N_t    = " << N_t << endl;
    cout << "N_step = " << N_step << endl;
    cout << "skip   = " << skip << endl;
    cout << "pausa  = " << pausa << endl;
    cout << "sigma  = " << sigma << endl;
    cout << "\n" << endl;
}
void stampa_coord(vec *r) {
    for (int i = 0; i < N; ++i) {
        cout << "#" << i << "\t" << r[i].x << "\t" << r[i].y << "\t" << r[i].z << "\t" << endl;
    }
    LOG("\n");
}
double pow1(double base, int esp) {//funzione esponenziale creata per non usare pow
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}

void plot_v_hist(vec *v) {
    //usare N grande
    double v_min = 1e100;
    double v_max = 0;
    double v_mod[N];
    for (int i = 0; i < N; ++i) {
        v_mod[i] = v[i].mod();
        if (v_mod[i] < v_min) {
            v_min = v_mod[i];
        } else if (v_mod[i] > v_max) {
            v_max = v_mod[i];
        }
    }
    int N_v = 1 + 3.322 * log(N); //Sturge's Rule
    double v_step = (v_max - v_min) / N_v;
    double bins;
    double f_v;
    for (int j = 0; j < N_v; ++j) {
        bins = v_min + v_step * j;
        f_v = 0;
        for (int i = 0; i < N; ++i) {
            if (v_mod[i] >= bins && v_mod[i] < (bins + v_step)) {
                f_v++;
            }
        }
        dati << bins << "\t" << f_v << endl;
    }
}
double min(double a, double b){
    if(a<=b){
        return a;
    }
    else{
        return b;
    }
}

#endif 