#ifndef es11_h
#define es11_h

#include "header.h"
#include "distribuzioni.h"

using namespace std;
using namespace arma;

extern int M;
extern int N;
extern int numero_accettati;
extern int numero_proposti;


double mod(cube r, double riga, double colonna){
    double mod= sqrt(pow1(r(riga,colonna,0),2)+pow1(r(riga,colonna,1),2)+pow1(r(riga,colonna,2),2));
    return mod;
}
double mod(mat r, double riga){
    double mod= sqrt(pow1(r(riga, 0),2)+pow1(r(riga, 1),2)+pow1(r(riga, 2),2));
    return mod;
}
void crea_reticolo(mat r, double L) {
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

double V_LJ(double r, double L) {
    //potenziale va a zero in modo continuo sul bordo della scatola
    if (r < L / 2 && r != 0) {
        double VL_2 = 4 * (pow1(2 / L, 12) - pow1(2 / L, 6));
        //double VpL_2 = 24 * (pow1(1 / r, 8) - 2 * pow1(1 / r, 14));
        return 4 * (pow1(1 / r, 12) - pow1(1 / r, 6)) - VL_2;
    } else {
        return 0;
    }
}

void accetto_spostamento(mat r, rowvec r_n, double V_tot_r0, double V_tot_r1, int n, double T){//accetto lo spostamento di MRT2?
    double A=min(1,exp(-(V_tot_r1-V_tot_r0)/T));//trovo A
    numero_proposti++;
    if(A>1){//se l'energia diminuisce accetto sempre lo spostamento
        r.row(n)=r_n;
        numero_accettati++;
    }
    else{//se l'energia aumenta accetto con probabilita uniforme come termostato di anderson
        if(A>(rand()/((double)RAND_MAX+1.0))){
            r.row(n)=r_n;
            numero_accettati++;
        }
        // else{//se non accetto lascio invariato
        //     r.row(n)=r.row(n);
        // }
    }
}
void posiz_MRT2(cube dr, cube dr_n, rowvec r_n, mat r, int i, double L, int n){
    for (int j = i + 1; j < N; ++j) {
        //calcolo la distanza tra la particella i e la particella j>i
        
        for(int k=0; k<3; ++k){//ciclo sulle coordinate
            dr(i,j,k) = r(i,k) - r(j,k);
            dr(j,i,k) = -dr(i,j,k);

            dr(i,j,k) -= L * rint(dr(i,j,k)/L);//sposto in [-L/2,+L/2]
            dr(j,i,k) -= L * rint(dr(j,i,k)/L);

            if(i!=n && j!=n){//calcolo le distanze nuove e impongo che siano le stesse se sono distanze tra particelle non modificate
                dr_n(i,j,k)=dr(i,j,k);
                dr_n(j,i,k)=-dr_n(i,j,k);
            }
            else if(i==n && j!=n){//modifico le distanze per una particella modificata
                dr_n(i,j,k) = r_n(k) - r(j,k);
                dr_n(j,i,k) = -dr_n(i,j,k);
                
                dr_n(i,j,k) -= L * rint(dr_n(i,j,k)/L);//sposto in [-L/2,+L/2]
                dr_n(j,i,k) -= L * rint(dr_n(j,i,k)/L);
            }
            else if (i!=n && j==n){//modifico le distanze per una particella modificata
                dr_n(i,j,k) = r(i,k) - r_n(k);
                dr_n(j,i,k) = -dr_n(i,j,k);
    
                dr_n(i,j,k) -= L * rint(dr_n(i,j,k)/L);//sposto in [-L/2,+L/2]
                dr_n(j,i,k) -= L * rint(dr_n(j,i,k)/L);
            }
            else{//impongo distanza tra particella modificata e se stessa uguale a zero
                dr_n(i,j,k) = 0;
            }
        }
    }
}

void MRT2(mat r, double *V, double *W, int N, double Delta, double T, double L, double r_c){//N=N_mol
    cube dr_n(N, N, 3), dr(N, N, 3);
    rowvec r_n; 
    
    int n = Fabs(rint(N * rand() / (RAND_MAX + 1.0)));//trovo la molecola che viene modificata da MTR2

    r_n(0) = r(n,0) + Delta * (rand() / (RAND_MAX + 1.0) - 0.5);//modifico le posizioni della particella n
    r_n(1) = r(n,1) + Delta * (rand() / (RAND_MAX + 1.0) - 0.5);
    r_n(2) = r(n,2) + Delta * (rand() / (RAND_MAX + 1.0) - 0.5);
    
    double V_tot_r1=0;//potenziale in posizione nuova
    double V_tot_r0=0;//potenziale in posizione vecchia

    double dr_mod, dr_mod_n;

    for (int i = 0; i < N; ++i) {

        posiz_MRT2(dr, dr_n, r_n, r, i, L, n);//sistema le posizioni vecchie e nuove in dr e dr_n

        for (int j = i + 1; j < N; ++j) {//trovo i potenziali nelle due posizioni
            
            dr_mod = mod(dr,i,j);
            dr_mod_n = mod(dr_n,i,j);
            
            if (dr_mod < r_c) {
                V_tot_r0+=V_LJ(dr_mod, L);
            }
            if (dr_mod_n < r_c) {
                V_tot_r1+=V_LJ(dr_mod_n, L);
            }
        }
    }
    
    accetto_spostamento(r, r_n, V_tot_r0, V_tot_r1, n, T);//verifico se accettare lo spostamento con MTR2
    
    *V = 0; *W = 0;
    for (int i = 0; i < N; ++i){//poco efficiente ma trovo e sparo fuori i valori di potenziale e pressione
        for (int j = i + 1; j < N; ++j) {
            dr_mod=mod(dr,i,j);
            if (dr_mod < r_c) {
                *V+=V_LJ(dr_mod, L);
                
                //dV_dr * r
                *W -= 24 * (pow1(1 / dr_mod, 6) - 2 * pow1(1 / dr_mod, 12));
            }
        }
    }
    *W/=N;
}
#endif 


