//
//  funzioni_utili.h
//  es1
//  Created by Mattia Lupi on 18/04/22.
//

#ifndef funzioni_utili_h
#define funzioni_utili_h
using namespace std;

extern ofstream dati;
extern int M;
double pow1(double base, int esp) {//funzione esponenziale creata per non usare pow
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}
long double pow1(long double base, int esp) {//funzione esponenziale creata per non usare pow versione long double
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}
void gauss_distr(double *v, double sigma, int N_Mol) { //crea distribuzione gaussiana per un vettore in ingresso
    //può essere ottimizzato: stiamo buttando via un valore per evitare che siano dispari
    for (int i = 0; i < 3 * N_Mol; i += 1) {
        double x1, x2;
        x1 = rand() / ((double)RAND_MAX + 1.0);
        x2 = rand() / ((double)RAND_MAX + 1.0);
        v[i] = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
        //v[i + 1] = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1);
    }
}

void setta_matr(double *r, double val, int N_Mol) {//set di un intero vettore ad un valore specificato
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = val;
    }
}
void setta_matr(long double *r, double val, int N_Mol) {//set di un intero vettore ad un valore specificato versione long double
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = val;
    }
}
void compila_matr(double *r, double *a_prev, int N_Mol,
                  void (*f)(double*, double, double*, int), double L) {//assegna alla matrice a_prev il valore di accelerazione imposto dal potenziale
    for (int j = 0; j < N_Mol; j++) {
        double F[3];
        f(r, L, F, j);//aggiorna F[3]
        for (int i = 3 * j; i < 3 * j + 3; ++i) {
            a_prev[i] = F[i - 3 * j];
        }
    }
}
void compila_matr(double *r, double *a_prev, int N_Mol,
                  void (*f)(double*, double, long double*, int), double L) { //assegna alla matrice a_prev il valore di accelerazione imposto dal potenziale versione long double
    for (int j = 0; j < N_Mol; j++) {
        long double F[3];
        f(r, L, F, j);//aggiorna F[3] usando fLJ
        for (int i = 3 * j + 0; i < 3 * j + 3; ++i) {
            a_prev[i] = F[i - 3 * j];
        }
    }
}
void copia_vett(double *r, int N_Mol, double *v) {
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = v[j];
    }
}
void set_vcm0(double *v, int N_Mol) {//imposta la velocita media a zero
    double v_cm[] = {0, 0, 0};
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < N_Mol; i++) {
            v_cm[j] += v[j + 3 * i];
        }
    }
    for (int j = 0; j < 3; j++) {
        v_cm[j] /= N_Mol;
    }
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < N_Mol; i++) {
            v[i + 3 * j] -= v_cm[j];
        }
    }
}
void crea_reticolo(int N_mol, double L, double *r) {
    if(M==1){
        double x, y, z, b;
        b = L / cbrt(N_mol);
        int k = 0; //contatore della particella
        for (x = 0; x < L; x += b) {
            for (y = 0; y < L; y += b) {
                for (z = 0; z < L; z += b) {
                    r[k] = x;
                    r[k + 1] = y;
                    r[k + 2] = z;
                    k+=3;
                }
            }
        }
    }
    else if(M==2){
        double a[]={0,0,0,0.5,0.5,0.5};
        double x, y, z, b;
        b = L / cbrt(N_mol/M);
        int k = 0; //contatore della particella
        for (x = 0; x < L; x += b) {
            for (y = 0; y < L; y += b) {
                for (z = 0; z < L; z += b) {
                    for(int j=0;j<6;j+=3){
                        r[k] = x+b*a[j];
                        r[k + 1] = y+b*a[j+1];
                        r[k + 2] = z+b*a[j+2];
                       k+=3;
                    }
                }
            }
        }
    }
    else if(M==4){
        double a[]={0,0,0,0.5,0.5,0,0.5,0,0.5,0,0.5,0.5};
        double x, y, z, b;
        b = L / cbrt(N_mol/M);
        int k = 0; //contatore della particella
        for (x = 0; x < L; x += b) {
            for (y = 0; y < L; y += b) {
                for (z = 0; z < L; z += b) {
                    for(int j=0;j<12;j+=3){
                        r[k] = x+b*a[j];
                        r[k + 1] = y+b*a[j+1];
                        r[k + 2] = z+b*a[j+2];
                        k+=3;
                    }
                }
            }
        }
    }
}
void stampa_reticolo(int N_mol, double *r) {
    for (int i = 0; i < 3 * N_mol; i += 3) {
        dati << r[i] << "\t" << r[i + 1] << "\t" << r[i + 2] << "\t" << i/3 << endl;
    }
    dati << "\n\n";
}
void cond_bordo(double *r, int N_mol, double L) { //impone condizioni di bordo periodiche su tutte le particelle
    for (int j = 0; j < 3 * N_mol; j++) {
        r[j] = r[j] - L * rint(r[j] / L);
    }
}


#endif /* funzioni_utili_h */