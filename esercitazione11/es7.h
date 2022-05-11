#include "header.h"
#include "distribuzioni.h"



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
        double VL_2 = 4 * (pow(2 / L, 12) - pow(2 / L, 6));
        //double VpL_2 = 24 * (pow(1 / r, 8) - 2 * pow(1 / r, 14));
        return 4 * (pow(1 / r, 12) - pow(1 / r, 6)) - VL_2;
    } else {
        return 0;
    }
}
void a_LJ(vec *r, vec *a, vec *dr[] , double r_c, double L) {
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            //calcolo la distanza tra la particella i e la particella j>i
            dr[i][j].x = r[i].x - r[j].x;
            dr[i][j].y = r[i].y - r[j].y;
            dr[i][j].z = r[i].z - r[j].z;

            dr[i][j].x -= L * rint(dr[i][j].x / L);
            dr[i][j].y -= L * rint(dr[i][j].y / L);
            dr[i][j].z -= L * rint(dr[i][j].z / L);
            dr[j][i].x -= L * rint(dr[j][i].x / L);
            dr[j][i].y -= L * rint(dr[j][i].y / L);
            dr[j][i].z -= L * rint(dr[j][i].z / L);
        }

        a[i].uguale(0);
    }
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dr_mod = dr[i][j].mod();
            if (dr_mod < r_c) {
                double cost = -24 * (pow(1 / dr_mod, 8) - 2 * pow(1 / dr_mod, 14));
                a[i].x += dr[i][j].x * cost;
                a[i].y += dr[i][j].y * cost;
                a[i].z += dr[i][j].z * cost;

                //sfrutto la forza uguale e contraria
                a[j].x += -dr[i][j].x * cost;
                a[j].y += -dr[i][j].y * cost;
                a[j].z += -dr[i][j].z * cost;

            }
        }
    }
}
void vel_verlet(vec *r, vec *v, vec *a, double dt, double r_c, double L, double *K, double *V, double *W) {
    vec a_prev[N], *dr[N]; vec_2D(dr, N);
    for (int i = 0; i < N; ++i) {
        //POSIZIONI
        r[i].x += v[i].x * dt + 0.5 * dt * dt * a[i].x;
        r[i].y += v[i].y * dt + 0.5 * dt * dt * a[i].y;
        r[i].z += v[i].z * dt + 0.5 * dt * dt * a[i].z;

        //sposto dentro alla scatola
        // r[i].x -= L * rint(r[i].x / L);
        // r[i].y -= L * rint(r[i].y / L);
        // r[i].z -= L * rint(r[i].z / L);

        // r[i].x += L/2;
        // r[i].y += L/2;
        // r[i].z += L/2;

        while (r[i].x > L) {r[i].x -= L;}
        while (r[i].y > L) {r[i].y -= L;}
        while (r[i].z > L) {r[i].z -= L;}

        while (r[i].x < 0) {r[i].x += L;}
        while (r[i].y < 0) {r[i].y += L;}
        while (r[i].z < 0) {r[i].z += L;}

        //copio le accelerazioni
        a_prev[i].x = a[i].x;
        a_prev[i].y = a[i].y;
        a_prev[i].z = a[i].z;
    }

    //ACCELERAZIONI
    a_LJ(r, a, dr, r_c, L);

    for (int i = 0; i < N; ++i) {
        //VELOCITA
        v[i].x += 0.5 * dt * (a_prev[i].x + a[i].x);
        v[i].y += 0.5 * dt * (a_prev[i].y + a[i].y);
        v[i].z += 0.5 * dt * (a_prev[i].z + a[i].z);
    }

    *K = 0; *V = 0; *W = 0;
    double v_mod, dr_mod;
    for (int i = 0; i < N; ++i) {
        v_mod = v[i].mod();
        *K += 0.5 * v_mod * v_mod;
        for (int j = i + 1; j < N; ++j) {
            dr_mod = dr[i][j].mod();
            if (dr_mod < r_c) {
                *V += V_LJ(dr_mod, L);

                //dV_dr * r
                *W -= 24 * (pow(1 / dr_mod, 6) - 2 * pow(1 / dr_mod, 12));
            }
        }
    }
    *W /= N;
}

void MRT2(vec *r, double *V, double *W, int N_mol, double Delta,){
    int n = rint(N_mol * rand() / (RAND_MAX + 1.0));//molecola che viene modificata da metropolis
    struct vec r_n, dr[N_mol], dr_n[N_mol];
    r_n.x=r[n].x+Delta*(rand() / (RAND_MAX + 1.0)-0.5);//modifico le posizioni della particella n
    r_n.y=r[n].y+Delta*(rand() / (RAND_MAX + 1.0)-0.5);
    r_n.z=r[n].z+Delta*(rand() / (RAND_MAX + 1.0)-0.5);

    V_tot_r1=0;
    V_tot_r0=0;
    double dr_mod;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            //calcolo la distanza tra la particella i e la particella j>i
            dr[i][j].x = r[i].x - r[j].x;
            dr[i][j].y = r[i].y - r[j].y;
            dr[i][j].z = r[i].z - r[j].z;
            dr[i][j].x -= L * rint(dr[i][j].x / L);
            dr[i][j].y -= L * rint(dr[i][j].y / L);
            dr[i][j].z -= L * rint(dr[i][j].z / L);
            dr[j][i].x -= L * rint(dr[j][i].x / L);
            dr[j][i].y -= L * rint(dr[j][i].y / L);
            dr[j][i].z -= L * rint(dr[j][i].z / L);

            if(i!=n && j!=n){
                dr_n[i][j].x=dr[i][j].x;
                dr_n[i][j].y=dr[i][j].y;
                dr_n[i][j].z=dr[i][j].z;
            }
        }
        for (int j = i + 1; j < N; ++j) {
            dr_mod = dr[i][j].mod();
            if (dr_mod < r_c) {
                V_tot_r0+=V_LJ(dr_mod, L);
                *V = V_tot_r0;

                //dV_dr * r
                *W -= 24 * (pow(1 / dr_mod, 6) - 2 * pow(1 / dr_mod, 12));
            }
        }
    }
    
    A=min(1,exp())


}















