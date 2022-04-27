//
//  funzioni_utili.h
//  es1
//  Created by Mattia Lupi on 18/04/22.
//

#ifndef funzioni_utili_h
#define funzioni_utili_h
using namespace std;

extern ofstream dati;

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
<<<<<<< Updated upstream
=======
<<<<<<< HEAD
void gauss_distr(double *v, double sigma, int N_Mol) {
    for (int i = 0; i < 3 * N_Mol; i += 2) {//tutte 3N coordinate
=======
>>>>>>> Stashed changes
void gauss_distr(double *v, double sigma, int N_Mol) { //crea distribuzione gaussiana per un vettore in ingresso
    for (int i = 0; i < 3 * N_Mol; i += 2) {
>>>>>>> a14900cb9df682e2da9c140b85653d8582dcf285
        double x1, x2;
        x1 = rand() / ((double)RAND_MAX + 1.0);
        x2 = rand() / ((double)RAND_MAX + 1.0);
        v[i] = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
        v[i + 1] = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1);
    }
}

void setta_matr(double *r, int val, int N_Mol) {//set di un intero vettore ad un valore specificato
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = val;
    }
}
void setta_matr(long double *r, int val, int N_Mol) {//set di un intero vettore ad un valore specificato versione long double
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = val;
    }
}
void compila_matr(double *r, double *a_prev, int N_Mol,
                  void (*f)(double*, double*, double*, int),
                  double *args) {//assegna alla matrice a_prev il valore di accelerazione imposto dal potenziale
    for (int j = 0; j < N_Mol; j++) {
        double F[3];
        f(r, args, F, j);//aggiorna F[3]
        for (int i = 3 * j; i < 3 * j + 3; ++i) {
            a_prev[i] = F[i - 3 * j];
        }
    }
}
void compila_matr(double *r, double *a_prev, int N_Mol,
<<<<<<< HEAD
                  void (*f)(double*, double*, long double*, int),
                  double *args) {
    for (int j = 0; j < N_Mol; j++) {
        long double F[3];
        f(r, args, F, j);//aggiorna F[3]
        for (int i = 3 * j; i < 3 * j + 3; ++i) {
            a_prev[i] = F[i - 3 * j];
        }
    }
}
=======
                   void (*f)(double*, double*, long double*, int),
                   double *args) {//assegna alla matrice a_prev il valore di accelerazione imposto dal potenziale versione long double
     for (int j = 0; j < N_Mol; j++) {
         long double F[3];
         f(r, args, F, j);//aggiorna F[3]
         for (int i = 3 * j; i < 3 * j + 3; ++i) {
             a_prev[i] = F[i - 3 * j];
         }
     }
 }
>>>>>>> a14900cb9df682e2da9c140b85653d8582dcf285
void copia_vett(double *r, int N_Mol, double *v) {
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = v[j];
    }
}
void set_vcm0(double *v, int N_Mol) {//imposta la velocita media a zero
    double v_cm[] = {0, 0, 0};
    for (int j = 0; j < 3; j++) {//per 3 direzioni
        for (int i = 0; i < N_Mol; i++) {//per N_mol particelle
            v_cm[j] += v[i + 3 * j];
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
    double x, y, z, b;
    b = L / cbrt(N_mol);
<<<<<<< Updated upstream
=======
<<<<<<< HEAD
    int k = 0; //contatore della posizione
    for (x = 0; x < L - b; x += b) {
        for (y = 0; y < L - b; y += b) {
            for (z = 0; z < L - b; z += b) {
                r[k] = x;
                r[k + 1] = y;
                r[k + 2] = z;
                k += 3;
                if (k > 3 * N_mol)
                    break;
=======
>>>>>>> Stashed changes
    int k = 0; //contatore della particella
    for (x = 0; x < L; x += b) {
        for (y = 0; y < L; y += b) {
            for (z = 0; z < L; z += b) {
                r[k] = x;
                r[k + 1] = y;
                r[k + 2] = z;
                //cout << r[k] << "\t" << r[k + 1] << "\t" << r[k + 2] << "\t" << k/3 << endl;
                dati << r[k] << "\t" << r[k + 1] << "\t" << r[k + 2] << "\t" << k/3 << endl;
                k+=3;
<<<<<<< Updated upstream
=======
>>>>>>> a14900cb9df682e2da9c140b85653d8582dcf285
>>>>>>> Stashed changes
            }
        }
    }
}
void stampa_reticolo(int N_mol, double *r) {
    for (int i = 0; i < 3 * N_mol; i += 3) {
        dati << r[i] << "\t" << r[i + 1] << "\t" << r[i + 2] << "\t" << i / 3<< endl;
    }
    dati << "\n\n";
}
<<<<<<< Updated upstream
=======
<<<<<<< HEAD
void cond_bordo(double *r, int N_mol, double L) {
    /*se il valore di una coordinata è maggiore di L,
      trasla la particella di -L*/
    for (int j = 0; j < 3 * N_mol; j++) {
        r[j] = r[j] - L * rint(r[j] / L);
=======
>>>>>>> Stashed changes
void cond_bordo(double *r,int N_mol, double L) {//impone condizioni di bordo periodiche su tutte le particelle
    for (int j = 0; j < 3*N_mol; j++) {
        r[j]=r[j]-L*rint(r[j]/L);
>>>>>>> a14900cb9df682e2da9c140b85653d8582dcf285
    }
}


#endif /* funzioni_utili_h */
