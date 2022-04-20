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
double fLJ(double *r, double *args, int indice, double r_mod);
void primi_vicini(double *r, double *r_prim_vic, double distanza_interaz, int N_mol, int j);

int main() {
    dati.open("dati.dat");
    double t=0, h, E, K, T;
    double t0 = 0, t1 = 25;
    int N_mol = pow1(10,3);
    int N = 1000;
    h = (t1 - t0) / ( (double) N);
    double rho=1e-2; //fisso la densita del campione da studiare
    double eps, sigma, L, distanza_interaz;
    L=cbrt(N_mol)/rho;
    
    double var_ad[]={eps, sigma, distanza_interaz, L};

    double r[3 * N_mol], v[3 * N_mol], r_prim_vic[3*N_mol];
    double r_mod, v_mod;
    
    r_mod = sqrt(r[0] * r[0] + r[0 + 1] * r[0 + 1] + r[0 + 2] * r[0 + 2]);
    v_mod = sqrt(v[0] * v[0] + v[0 + 1] * v[0 + 1] + v[0 + 2] * v[0 + 2]);
    
    crea_reticolo(N_mol, L, r);

    gauss_distr(v, 1, N_mol);
    set_vcm0(v,N_mol);

    double a_prev[3 * N_mol];
    compila_matr(r, a_prev, N_mol, fLJ, var_ad);//inizializzo accelerazioni come date dal potenziale

    E = 0; K = 0;
    for (int n = 0; n < N; ++n) {//passi temporali
        //E = E*(n+1);//energia totale @t
        K = K*(n+1);//energia cinetica @t
        T = 0;//temp @t
        
        r_mod = sqrt(r[0] * r[0] + r[0 + 1] * r[0 + 1] + r[0 + 2] * r[0 + 2]);
        v_mod = sqrt(v[0] * v[0] + v[0 + 1] * v[0 + 1] + v[0 + 2] * v[0 + 2]);
        
        for (int j = 0; j < N_mol; ++j) {//particelle
            t = t0 + n * h;
            
            primi_vicini(r,r_prim_vic,var_ad[2],N_mol,j);
            if (vel_verlet(t, r, v, h, j, 3, a_prev, fLJ, var_ad, r_mod)) {printf("ERRORE");}

            r_mod = sqrt(r[j] * r[j] + r[j + 1] * r[j + 1] + r[j + 2] * r[j + 2]);
            v_mod = sqrt(v[j] * v[j] + v[j + 1] * v[j + 1] + v[j + 2] * v[j + 2]);
            
            cond_bordo(r,j,L);
            
            
            K += 0.5 * v_mod * v_mod;
            //E += 0.5 * r_mod * r_mod + 0.5 * v_mod * v_mod;
            
        }
       // E = E/(n+2);
        K = K/(n+2);
        T = 2.0*K/(3.0*N_mol);
        //dati << t  << "\t" << E << "\t" << K << "\t" << T << endl;
    }
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
void fLJ(double *r, double *args, double *f, int N_mol, int i){//arg[0]=eps, arg[1]=sigma, arg[2]= dimensione scatola per particella, arg[3]=dimensione scatola
    double r_pv[3*N_mol];
    primi_vicini(r,r_pv,args[3],N_mol,i);
    double modR_Prim_vic[N_mol];
    for(int j=0; j<3*N_mol; j+=3){
        modR_Prim_vic[j/3]=sqrt(r_pv[j]*r_pv[j]+r_pv[j+1]*r_pv[j+1]+r_pv[j+2]*r_pv[j+2]);
    }
    for(int j=0; j<3*N_mol; j++){
        if((r_pv[j])>=args[2]/2){//manda la forza a zero se piu distante di dist/2
            switch (j%3) {
                case 0:
                    f[0]+=0;
                    break;
                case 1:
                    f[1]+=0;
                    break;
                case 2:
                    f[2]+=0;
                    break;
            }
        }
    }
    else{
        if(j!=i && j!=i+1 && j!=i+2){
            double r7=pow1(modR_Prim_vic[i],7);
            double sigma6=pow1(args[1],6);
            f[=24*args[0]*(2*r_pv[j]*sigma6*sigma6/(r7*r7)-sigma6*r[i]/(r7*r_mod));
        }
    }
}
void primi_vicini(double *r, double *r_prim_vic, double distanza_interaz, int N_mol, int j){
    for(int i=0; i<3*N_mol; i+=3){
        r_prim_vic[i]=r[i]+distanza_interaz*round((r[j]-r[i])/distanza_interaz);
        r_prim_vic[i+1]=r[i+1]+distanza_interaz*round((r[j+1]-r[i+1])/distanza_interaz);
        r_prim_vic[i+2]=r[i+2]+distanza_interaz*round((r[j+2]-r[i+2])/distanza_interaz);
    }
}

