#ifndef es11_h
#define es11_h

#include "header.h"
#include "distribuzioni.h"

using namespace std;

extern int M;
extern int N;

extern int numero_accettati;
extern int numero_proposti;

void crea_reticolo(vec *r, double L) {
    int n = cbrt(N / M);
    float L_cella = L / cbrt(N / M);
    int cont = 0;//contatore particella

    double b[][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    if (M == 2) {b[1][2] = 0.5;}
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                vec R = {.x = i * L_cella, .y = j * L_cella, .z = k * L_cella};
                for (int l = 0; l < M; ++l) {
                    r[cont].x = R.x + b[l][0] * L_cella;
                    r[cont].y = R.y + b[l][1] * L_cella;
                    r[cont].z = R.z + b[l][2] * L_cella;
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

void accetto_spostamento(vec *r, vec r_n, double V_tot_r0, double V_tot_r1, int n, double T){//accetto lo spostamento di MRT2?
    double A=min(1,exp(-(V_tot_r1-V_tot_r0)/T));//trovo A
    numero_proposti++;
    if(A>1){//se l'energia diminuisce accetto sempre lo spostamento
        r[n].x=r_n.x;
        r[n].y=r_n.y;
        r[n].z=r_n.z;
        numero_accettati++;
    }
    else{//se l'energia aumenta accetto con probabilita uniforme come termostato di anderson
        if(A>(rand()/((double)RAND_MAX+1.0))){
            r[n].x=r_n.x;
            r[n].y=r_n.y;
            r[n].z=r_n.z;
            numero_accettati++;
        }
        // else{
        //     r[n].x=r[n].x;
        //     r[n].y=r[n].y;
        //     r[n].z=r[n].z;
        // }
    }
}
void posiz_MRT2(vec *dr[], vec *dr_n[], vec r_n, vec *r, int i, double L, int n){
    for (int j = i + 1; j < N; ++j) {
        //calcolo la distanza tra la particella i e la particella j>i
        dr[i][j].x = r[i].x - r[j].x;
        dr[i][j].y = r[i].y - r[j].y;
        dr[i][j].z = r[i].z - r[j].z;
        dr[i][j].x -= L * rint(dr[i][j].x / L);//sposto in [-L/2,+L/2]
        dr[i][j].y -= L * rint(dr[i][j].y / L);
        dr[i][j].z -= L * rint(dr[i][j].z / L);
        dr[j][i].x -= L * rint(dr[j][i].x / L);
        dr[j][i].y -= L * rint(dr[j][i].y / L);
        dr[j][i].z -= L * rint(dr[j][i].z / L);
        if(i!=n && j!=n){//calcolo le distanze nuove e impongo che siano le stesse se sono distanze tra particelle non modificate
            dr_n[i][j].x=dr[i][j].x;
            dr_n[i][j].y=dr[i][j].y;
            dr_n[i][j].z=dr[i][j].z;
        }
        else if(i==n && j!=n){//modifico le distanze per una particella modificata
            dr_n[i][j].x = r_n.x - r[j].x;
            dr_n[i][j].y = r_n.y - r[j].y;
            dr_n[i][j].z = r_n.z - r[j].z;
            dr_n[i][j].x -= L * rint(dr_n[i][j].x / L);//sposto in [-L/2,+L/2]
            dr_n[i][j].y -= L * rint(dr_n[i][j].y / L);
            dr_n[i][j].z -= L * rint(dr_n[i][j].z / L);
            dr_n[j][i].x -= L * rint(dr_n[j][i].x / L);
            dr_n[j][i].y -= L * rint(dr_n[j][i].y / L);
            dr_n[j][i].z -= L * rint(dr_n[j][i].z / L);
        }
        else if (i!=n && j==n){//modifico le distanze per una particella modificata
            dr_n[i][j].x = r[i].x - r_n.x;
            dr_n[i][j].y = r[i].y - r_n.y;
            dr_n[i][j].z = r[i].z - r_n.z;
            dr_n[i][j].x -= L * rint(dr_n[i][j].x / L);//sposto in [-L/2,+L/2]
            dr_n[i][j].y -= L * rint(dr_n[i][j].y / L);
            dr_n[i][j].z -= L * rint(dr_n[i][j].z / L);
            dr_n[j][i].x -= L * rint(dr_n[j][i].x / L);
            dr_n[j][i].y -= L * rint(dr_n[j][i].y / L);
            dr_n[j][i].z -= L * rint(dr_n[j][i].z / L);
        }
        else{//impongo distanza tra particella modificata e se stessa uguale a zero
            dr_n[i][j].uguale(0);
        }
    }
}

void MRT2(vec *r, double *V, double *W, int N, double Delta, double T, double L, double r_c){//N=N_mol
    vec *dr_n[N], *dr[N], r_n; 
    vec_2D(dr, N);
    vec_2D(dr_n,N);//creo le matrici di struct

    int n = Fabs(rint(N * rand() / (RAND_MAX + 1.0)));//trovo la molecola che viene modificata da MTR2

    r_n.x=r[n].x+Delta*(rand() / (RAND_MAX + 1.0)-0.5);//modifico le posizioni della particella n
    r_n.y=r[n].y+Delta*(rand() / (RAND_MAX + 1.0)-0.5);
    r_n.z=r[n].z+Delta*(rand() / (RAND_MAX + 1.0)-0.5);
    
    double V_tot_r1=0;//potenziale in posizione nuova
    double V_tot_r0=0;//potenziale in posizione vecchia

    double dr_mod, dr_mod_n;

    for (int i = 0; i < N; ++i) {

        posiz_MRT2(dr, dr_n, r_n, r, i, L, n);//sistema le posizioni vecchie e nuove in dr e dr_n

        for (int j = i + 1; j < N; ++j) {//trovo i potenziali nelle due posizioni
            
            dr_mod = dr[i][j].mod();
            dr_mod_n = dr_n[i][j].mod();
            
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
            dr_mod=dr[i][j].mod();
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


