#ifndef funzioni_utili_h
#define funzioni_utili_h

using namespace std;
extern ofstream dati;

void set_vcm0(double v[][3], int N_Mol){
    double v_cm[]={0,0,0};
    for(int j=0;j<3;j++){
        for(int i=0;i<N_Mol;i++){
            v_cm[j]+=v[i][j];
        }
    }
    for(int j=0;j<3;j++){
        v_cm[j]/=N_Mol;
    }
    for(int j=0;j<3;j++){
        for(int i=0;i<N_Mol;i++){
            v[i][j]-=v_cm[j];
        }
    }
}

void verifica_gauss(double *v, int N_bin, int N_Mol){
    for (double j = -20; j < 20; j += 40.0 / N_bin) {
        int cont = 0;
        for (int i = 0; i < N_Mol; i++) {
            if (v[i] < j + 40.0 / N_bin && v[i] > j) {
                cont++;
            }
        }
        dati << j << "\t" << cont << endl;
    }
}
void gauss_distr(double v[][3], double sigma, int N_Mol){
    for(int j=0;j<3;j++){
        for (int i = 0; i < N_Mol; i += 2) {
            double x1, x2;
            x1 = rand() / ((double)RAND_MAX + 1.0);
            x2 = rand() / ((double)RAND_MAX + 1.0);
            v[i][j] = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
            v[i + 1][j] = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1);
            //dati << i << "\t" << v[i] << endl;
        }
    }
}
void resetta_matr(double r[][3], int N_Mol){
    for(int j=0;j<3;j++){
        for (int i = 0; i < N_Mol; i++) {
            r[i][j] = 0;
        }
    }
}

#endif /* funzioni_utili_h */
