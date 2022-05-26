#ifndef funz_gdr_h
#define funz_gdr_h

using namespace std;
using namespace arma;

extern int N;
extern int M;


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
void crea_reticolo_unif(mat &r, double L) {// passo la matrice per riferimento 
    double L_cella = L / cbrt(N / M);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 3; ++j) {
            r(i,j)=L*rand()/((double)RAND_MAX+1.0);
        }
        //cout<<i<<endl;
    }
}
void gdr_funz(mat &r, double L, double rho, mat &gdr, int N_b){
    double delta_r=L/((double)N_b*2);
    rowvec dr(3);
    double dr_mod;

    for (int k = 0; k < N_b; ++k){//ciclo sui bin
        for (int i = 0; i < N; ++i){//ciclo sulle particelle centrali
            double R=delta_r*k+delta_r/2;//definisco il raggio medio del volumetto sferico
            double dV=4/3*M_PI*(pow1(R+delta_r/2,3)-pow1(R-delta_r/2,3));//definisco il volumetto sferico
            double freq=0;//numero di particelle in un volumetto

            for (int j = i+1; j < N; ++j){//ciclo sulle particelle non centrali
                for (int k1 = 0; k1 < 3; ++k1){
                    dr(k1) = r(i,k1) - r(j,k1);
                }
                dr_mod=mod(dr);
                if(dr_mod<=R+delta_r/2 && dr_mod>R-delta_r/2){//se nel volumetto allora aumento la freq 
                    freq++;
                }
            }
            gdr(k,0)+=1/rho*(freq/dV);//def gdr
            gdr(k,1)=R;
        }
    }
    gdr.col(0)/=N;//medio sulle part
}


#endif