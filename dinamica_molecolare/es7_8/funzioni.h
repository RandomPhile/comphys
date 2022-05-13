#ifndef FUNZIONI_H
#define FUNZIONI_H

#include <iostream>
#include <fstream>
#include <cmath>
#define _USE_MATH_DEFINES
#define LOG(x) cout<<x<<endl;
using namespace std;

/*** numero di particelle e tipo di reticolo ***/
extern const int N;
extern const int M;
extern const double T_req;
struct coppia {
    double sigma; double rho;
};
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
    void stampa() {
        cout << "x: " << x << "y: " << y << "z: " << z << "\n";
    }
};

void crea_reticolo(vec *r, double L);

void distr_gauss(vec *x, double sigma);

int input();

void vec_2D (vec *dr[], int dim);

void aggiorna_a(vec *r, vec *a, double L);

void stampa_coord(vec *r, vec *v, ofstream &file);

void v_cm_0(vec *v);

void vel_verlet(vec *r, vec *v, vec *a, double dt, double L);

void calcolo_coordinate(string file_path, double rho, double sigma, double dt, double t1);

void plot_coordinate(string file_path, int N_step, double pausa);

void calcolo_osservabili_da_file(string coord_path, string obs_path);

void plot_osservabili();

double pressione(double rho, double sigma, double dt, double t1);

void calcolo_pressioni(string p_path, coppia *coppie, double dt, double t1, int c_min, int c_max);

void plot_pressioni();

void calcolo_coordinate_per_gdr(string coord_g_path, double rho, double sigma, double dt, double t1);

void calcolo_gdr_da_file(string coord_g_path, string g_path, double rho, int N_bins);

void plot_gdr();

#endif