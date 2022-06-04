#ifndef FUNZIONI_H
#define FUNZIONI_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include <armadillo>
#include <iomanip>
#include <algorithm>

// #define _USE_MATH_DEFINES
#define LOG(x) cout<<x<<endl;
using namespace std;
using namespace arma;


struct coppia {
    double sigma; double rho; double t_eq;
};

double pow1(double base, int esp) {//funzione esponenziale creata per non usare pow
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}

double mod(cube r, double riga, double colonna){//calcolo il modulo della posizione relativa delle particelle i e j
    double mod= sqrt(pow1(r(riga,colonna,0),2)+pow1(r(riga,colonna,1),2)+pow1(r(riga,colonna,2),2));
    return mod;
}
double mod(rowvec r){//calcolo il modulo della posizione relativa delle particelle i e j
    double mod= sqrt(pow1(r(0),2)+pow1(r(1),2)+pow1(r(2),2));
    return mod;
}
void crea_reticolo(mat &r, double L) {// passo la matrice per riferimento 
    int n = cbrt(N / M);
    double L_cella = L / cbrt(N / M);
    int cont = 0;//contatore particella

    mat b = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    if (M == 2) {b(1,2) = 0.5;}
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                rowvec R = {L_cella*(double)i, L_cella*(double)j, L_cella*(double)k};

                for (int l = 0; l < M; ++l) {
                    r.row(cont) = R + b.row(l) * L_cella;
                    cont++;
                }
            }
        }
    }
}
void distr_gauss(mat &x, double sigma, int num_righe, int num_col) {
    //Formule di Box-Muller
    for (int i = 0; i < num_righe; ++i) {
        for (int j = 0; j < num_col; ++j){
            double x1, x2;
            x1 = rand() / ((double)RAND_MAX + 1.0);
            x2 = rand() / ((double)RAND_MAX + 1.0);
    
            x(i,j) = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
            if(j!=num_col-1){
                x(i,j+1) = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1);
            }
        }
    }
}
int input() {
    int input;
    cout << "\n"
         "1. calcola coordinate per un valore di densità\n"
         "2. plotta coordinate da file\n"
         "3. calcola osservabili per un valore di densità\n"
         "4. plotta osservabili da file\n";
    cout << "step: "; cin >> input;
    return input;
}
void aggiorna_a(mat &r, mat &a, double L) {
    rowvec dr(3);
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < 3; ++k){
            a(i,k)=0;
        }

        for (int j = 0; j < N; ++j) {
            if (j != i) {
                for (int k = 0; k < 3; ++k){
                    dr(k) = r(i,k) - r(j,k);
                
                    dr(k) -= L * rint(dr(k) / L);
                }

                double dr_mod = mod(dr);

                if (dr_mod < L / 2) {
                    double cost = -24 * (pow(1 / dr_mod, 8) - 2 * pow(1 / dr_mod, 14));

                    for (int k = 0; k < 3; ++k){
                        a(i,k) += dr(k) * cost;
                    }
                }
            }
        }
    }
}
void stampa_coord(mat &r, mat &v, ofstream &file) {
    /* struttura file .xyz per VMD
    N
    nome molecola
    atomo1 x y z
    atomo2 x y z
    ...
    atomoN x y z
    */
    for (int i = 0; i < N; ++i) {
        file << "P" << i << "\t" << r(i,0) << "\t" << r(i,1) << "\t" << r(i,2) << "\t" <<v(i,0) << "\t" << v(i,1) << "\t" << v(i,2) << "\n";
    }
}
void v_cm_0(mat &v) {
    rowvec v_cm = {0, 0, 0};
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 3; ++j){
            v_cm(j) += v(i,j);
        }
    }
    v_cm/=N;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 3; ++j){
            v(i,j) -= v_cm(j);
        }
    }
}
void vel_verlet(mat &r, mat &v, mat &a, double dt, double L) {
    mat a_prev(N,3);
    for (int i = 0; i < N; ++i) {
        //POSIZIONI
        for (int j = 0; j < 3; ++j){
            r(i,j)+= v(i,j) * dt + 0.5 * dt * dt * a(i,j);

            while (r(i,j) > L) {r(i,j) -= L;}//sposto dentro alla scatola
            while (r(i,j) < 0) {r(i,j) += L;}

            //copio le accelerazioni
            a_prev(i,j) = a(i,j);
        }
    }

    //ACCELERAZIONI
    aggiorna_a(r, a, L);

    //VELOCITA
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 3; ++j){
            v(i,j) += 0.5 * dt * (a_prev(i,j) + a(i,j));
        }
    }
}
void calcolo_coordinate(string coord_path, double rho, double sigma, double dt, double t1) {
    double L = cbrt(N / rho);
    const int N_t = t1 / dt;
    mat r(N,3), v(N,3), a(N,3);//, *dr[N]; vec_2D(dr, N);
    ofstream coord;

    /*** inizializzo variabili ***/
    crea_reticolo(r, L);
    distr_gauss(v, sigma, N, 3);
    aggiorna_a(r, a, L);
    

    coord.open(coord_path);
    coord << N << "\t" << N_t << "\t" << dt << "\t" << rho << "\t" << sigma << "\n";
    coord << "molecola\n";

    for (int i = 0; i < N_t; ++i) {
        stampa_coord(r, v, coord);
        vel_verlet(r, v, a, dt, L);
    }

    coord.close();
}
void plot_coordinate(string coord_path, int N_step, double pausa) {
    ifstream coord_in; string comando;
    int N_t; double rho, sigma, dt;

    coord_in.open(coord_path);
    coord_in >> N_t >> N_t >> dt >> rho >> sigma;
    coord_in.close();

    if (N_step > N_t) {N_step = N_t;}
    int skip = rint(N_t / N_step); if (skip == 0) {skip = 1;}
    double L = cbrt(N / rho);

    comando = "gnuplot -e N=";
    comando += to_string(N);
    comando += " -e L=";
    comando += to_string(L);
    comando += " -e N_step=";
    comando += to_string(N_step);
    comando += " -e dt=";
    comando += to_string(dt);
    comando += " -e skip=";
    comando += to_string(skip);
    comando += " -e pausa=";
    comando += to_string(pausa);
    comando += " animation.plt";
    LOG(comando);
    system(comando.c_str());
}
void calcolo_osservabili_da_file(string coord_path, string obs_path, double t_eq) {
    ifstream coord_in; ofstream obs;
    string line;
    int N_t; double rho, sigma, dt;

    coord_in.open(coord_path);
    coord_in >> N_t >> N_t >> dt >> rho >> sigma;
    coord_in >> line; //scarta la seconda riga

    double L = cbrt(N / rho);

    mat r(N,3); 
    rowvec v(3), dr(3);
    double dr_mod, v_mod;

    obs.open(obs_path);
    double t = 0;
    double K_avg = 0, W_avg = 0, P, T;
    double K_avg_r = 0, T_r;
    double K, V, E, W;
    for (int i_t = 0; i_t < N_t; ++i_t) {

        K = 0; V = 0; W = 0;

        for (int i = 0; i < N; ++i) {
            coord_in >> line; //scarta la prima colonna

            coord_in >> r(i,0) >> r(i,1) >> r(i,2);
            coord_in >> v(0) >> v(1) >> v(2);

            v_mod = mod(v);
            
            K += 0.5 * v_mod * v_mod;
        }
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                for (int k = 0; k < 3; ++k){
                    dr(k) = r(i,k) - r(j,k);
                    dr(k) -= L * rint(dr(k)/L);//sposto in [-L/2,+L/2]
                }

                dr_mod = mod(dr);
                if (dr_mod < L / 2) {
                    if (dr_mod < L / 2 && dr_mod != 0) {
                        double VL_2 = 4 * (pow(2 / L, 12) - pow(2 / L, 6));
                        //double VpL_2 = 24 * (pow(1 / dr_mod, 8) - 2 * pow(1 / dr_mod, 14));
                        V += 4 * (pow(1 / dr_mod, 12) - pow(1 / dr_mod, 6)) - VL_2;
                    }
                    W += -24 * (pow(1 / dr_mod, 6) - 2 * pow(1 / dr_mod, 12));
                }
            }
        }
        W /= N;
        if(t>t_eq){
            K_avg_r = K_avg_r + (K - K_avg_r) / ((i_t - (int)(t_eq / dt)) + 1);//forse è +2
            T_r = 2.0 * K_avg_r / (3.0 * N);
        }


        K_avg = K_avg + (K - K_avg) / (i_t + 1);//forse è +2
        W_avg = W_avg + (W - W_avg) / (i_t + 1);//forse è +2

        T = 2.0 * K_avg / (3.0 * N);
        P = (1 + W_avg / (3.0 * T_req));

        E = K + V;
        obs << i_t << "\t" << K << "\t" << V << "\t" << E << "\t" << T << "\t" << P << "\n";
    
        t += dt;
    }
    cout << "Rho = " << rho << "\t\tT = " << T_r << "\nSigma = " << sigma << "\t->\t" << sigma*sqrt(T_req / T_r) << "\n" << endl;
    coord_in.close();
    obs.close();
}
void plot_osservabili() {
    string comando;

    comando = "gnuplot";
    comando += " plot.plt";
    LOG(comando);
    system(comando.c_str());
}
double pressione(double rho, double sigma, double dt, double t1, double t_eq) {
    double L = cbrt(N / rho);
    const int N_t = t1 / dt;
    mat r(N,3), v(N,3), a(N,3);//, *dr[N]; vec_2D(dr, N);
    rowvec dr(3);
    double dr_mod;

    /*** inizializzo variabili ***/
    crea_reticolo(r, L);
    distr_gauss(v, sigma, N, 3);
    aggiorna_a(r, a, L);

    double W_avg = 0, P;
    double W_avg_r = 0, P_r;
    double W;

    double t = 0;
    for (int i_t = 0; i_t < N_t; ++i_t) {
        W = 0;
        vel_verlet(r, v, a, dt, L);

        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                for (int k = 0; k < 3; ++k){
                    dr(k) = r(i,k) - r(j,k);
                    dr(k) -= L * rint(dr(k)/L);//sposto in [-L/2,+L/2]
                }

                dr_mod = mod(dr);
                if (dr_mod < L / 2) {
                    W += -24 * (pow(1 / dr_mod, 6) - 2 * pow(1 / dr_mod, 12));
                }
            }
        }
        W /= N;
        if(i_t >= (int)(t_eq / dt)){
            W_avg = W_avg + (W - W_avg) / (i_t + 1);//forse è +2
        }
    }
    return (1 + W_avg / (3.0 * T_req));
}
void calcolo_pressioni(string p_path, coppia *coppie, double dt, double t1, int c_min, int c_max, double t_eq) {
    ofstream press;
    press.open(p_path);
    for (int caso = c_min; caso <= c_max; ++caso) {
        press << coppie[caso].rho << "\t" << pressione(coppie[caso].rho, coppie[caso].sigma, dt, t1, t_eq) << "\n";
    }
    press.close();
}
void plot_pressioni() {
    string comando;

    comando = "gnuplot";
    comando += " plot_pressione.plt";
    LOG(comando);
    system(comando.c_str());
}
void calcolo_coordinate_per_gdr(string coord_g_path, double rho, double sigma, double dt, double t1) {
    ofstream coord_gdr;
    coord_gdr.open(coord_g_path);

    double L = cbrt(N / rho);
    const int N_t = t1 / dt;
    mat r(N,3), v(N,3), a(N,3);
    rowvec rp(3);
    int p = 0;//particella attorno cui calcolo g(r)

    /*** inizializzo variabili ***/
    crea_reticolo(r, L);
    distr_gauss(v, sigma, N, 3);
    aggiorna_a(r, a, L);

    double t = 0;
    for (int i_t = 0; i_t < N_t; ++i_t) {
        vel_verlet(r, v, a, dt, L);
    }
    for (int p = 0; p < N; ++p){
        coord_gdr << "Particella" << p << endl;
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k){
                rp(k) = r(p,k) - r(i,k);//trovo distanza relativa part p,j
    
                rp(k) -= L * rint(rp(k) / L);//sposto in [-L/2,+L/2]
            }
            coord_gdr << mod(rp) << "\n";
        }
    }
    coord_gdr.close();
}
void calcolo_gdr_da_file(string coord_g_path, string g_path, double rho, int N_bins) {
    ifstream coord_gdr;
    string line;
    coord_gdr.open(coord_g_path);
    ofstream gdr;

    double L = cbrt(N / rho);

    double delta_r = 0.5 * L / N_bins;
    cout << "delta_r = " << delta_r << endl;
    cout << "L/2 = " << L / 2 << endl;

    int freq[N_bins];
    for (int k = 0; k < N_bins; ++k) {
        freq[k] = 0;
    }
    double g, rp_mod, r_k;
    for (int p = 0; p < N; ++p){
        coord_gdr >> line;
        // cout<<line<<endl;
        for (int i = 0; i < N; ++i) {
            coord_gdr >> rp_mod;
            // cout<<rp_mod<<endl;
            for (int k = 0; k < N_bins; ++k) {
                r_k = (2 * k + 1) * delta_r / 2;
                if (rp_mod > r_k - delta_r / 2 && rp_mod <= r_k + delta_r / 2) {
                    freq[k]++;
                    break;
                }
            }
        }
    }
    coord_gdr.close();
    gdr.open(g_path);
    for (int k = 0; k < N_bins; ++k) {
        r_k = (2 * k + 1) * delta_r / 2;
        double dV = 4 * M_PI * r_k * r_k * delta_r + M_PI / 3 * pow1(delta_r, 3);//definisco il volumetto sferico
        g = freq[k] / dV;
        g /= rho;
        g /= N;
        gdr << r_k << "\t" << freq[k] << "\t" << g << "\n";
    }
    gdr.close();
}
void plot_gdr() {
    string comando;

    comando = "gnuplot";
    comando += " plot_gdr.plt";
    LOG(comando);
    system(comando.c_str());
}

 #endif