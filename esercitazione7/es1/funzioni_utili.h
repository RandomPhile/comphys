//
//  funzioni_utili.h
//  es1
//  Created by Mattia Lupi on 18/04/22.
//

#ifndef funzioni_utili_h
#define funzioni_utili_h
using namespace std;

extern ofstream dati;

double pow1(double base, int esp) {
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}
long double pow1(long double base, int esp) {
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}
void gauss_distr(double *v, double sigma, int N_Mol) {
    for (int i = 0; i < 3 * N_Mol; i += 2) {
        double x1, x2;
        x1 = rand() / ((double)RAND_MAX + 1.0);
        x2 = rand() / ((double)RAND_MAX + 1.0);
        v[i] = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
        v[i + 1] = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1);
    }
}

void setta_matr(double *r, int val, int N_Mol) {
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = val;
    }
}
void setta_matr(long double *r, int val, int N_Mol) {
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = val;
    }
}
void compila_matr(double *r, double *a_prev, int N_Mol,
                  void (*f)(double*, double*, double*, int),
                  double *args) {
    for (int j = 0; j < N_Mol; j++) {
        double F[3];
        f(r, args, F, j);//aggiorna F[3]
        for (int i = 3 * j; i < 3 * j + 3; ++i) {
            a_prev[i] = F[i - 3 * j];
        }
    }
}
void compila_matr(double *r, double *a_prev, int N_Mol,
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
void copia_vett(double *r, int N_Mol, double *v) {
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = v[j];
    }
}
void set_vcm0(double *v, int N_Mol) {
    double v_cm[] = {0, 0, 0};
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < N_Mol; i++) {
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
    int k = 0; //contatore della posizione
    for (x = 0; x < L-b; x += b) {
        for (y = 0; y < L-b; y += b) {
            for (z = 0; z < L-b; z += b) {
                r[k] = x;
                r[k + 1] = y;
                r[k + 2] = z;
                //cout << r[k] << "\t" << r[k + 1] << "\t" << r[k + 2] << "\t" << k/3 << endl;
                dati << r[k] << "\t" << r[k + 1] << "\t" << r[k + 2] << "\t" << k/3 << endl;
                k+=3;
                if (k > 3 * N_mol)
                    break;
            }
        }
    }
    dati<<"\n\n";
}
void stampa_reticolo(int N_mol, double L, double *r) {
    double x, y, z, b;
    b = L / cbrt(N_mol);
    int k = 0; //contatore della posizione
    for (x = 0; x < L; x += b) {
        for (y = 0; y < L; y += b) {
            for (z = 0; z < L; z += b) {
                dati << r[k] << "\t" << r[k + 1] << "\t" << r[k + 2] << "\t" << k << endl;
                //cout << r[k] << "\t" << r[k + 1] << "\t" << r[k + 2] << "\t" << k << endl;
                k+=3;
                if (k > 3 * N_mol)
                    break;
            }
        }
    }
}
void cond_bordo(double *r, int i, double L) {
    for (int j = i; j < i + 3; j++) {
        r[j]=r[j]-L*rint(r[j]/L);
    }
}


#endif /* funzioni_utili_h */
