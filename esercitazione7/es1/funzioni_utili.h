//
//  funzioni_utili.h
//  es1
//  Created by Mattia Lupi on 18/04/22.
//

#ifndef funzioni_utili_h
#define funzioni_utili_h

using namespace std;

extern ofstream dati;

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
void compila_matr(double *r,double *a, int N_Mol,
                  double (*f)(double*, double*, int, double),
                  double *fargs) {
    double r_mod;
    for (int j = 0; j < 3 * N_Mol; j++) {
        if(j%3==0){
            r_mod=sqrt(r[j] * r[j] + r[j + 1] * r[j + 1] + r[j + 2] * r[j + 2]);
        }
        else if(j%3==1){
            r_mod=sqrt(r[j] * r[j] + r[j + 1] * r[j + 1] + r[j - 1] * r[j - 1]);
        }
        else{
            r_mod=sqrt(r[j] * r[j] + r[j - 1] * r[j - 1] + r[j - 2] * r[j - 2]);
        }
        a[j] = f(r,fargs,j, r_mod);
    }
}
void copia_vett(double *r, int N_Mol,double *v) {
    for (int j = 0; j < 3 * N_Mol; j++) {
        r[j] = v[j];
    }
}
void set_vcm0(double *v, int N_Mol){
    double v_cm[]={0,0,0};
    for(int j=0;j<3;j++){
        for(int i=0;i<N_Mol;i++){
            v_cm[j]+=v[i+3*j];
        }
    }
    for(int j=0;j<3;j++){
        v_cm[j]/=N_Mol;
    }
    for(int j=0;j<3;j++){
        for(int i=0;i<N_Mol;i++){
            v[i+3*j]-=v_cm[j];
        }
    }
}
void crea_reticolo(int N_mol, double L, double *r){
    double x,y,z,b;
    b=L/cbrt(N_mol);
    int k=0;//contatore della posizione
    for(x=0; x<L; x+=b){
        for(y=0; y<L; y+=b){
            for(z=0; z<L; z+=b){
                r[k]=x;
                r[k+1]=y;
                r[k+2]=z;
                dati<<r[k]<<"\t"<<r[k+1]<<"\t"<<r[k+2]<<"\t"<<k<<endl;
                k++;
                if(k>3*N_mol)
                    break;
            }
        }
    }
}
void cond_bordo(double *r, int i, double L){
    for(int j=i; j<i+3; j++){//occhio a effetti brutti di spostamenti maggiori di L in tempo dt
        if(r[j]>L){
            r[j]-=L;
        }
        else if(r[j]<0){
            r[j]+=L;
        }
        else{}
    }
}
void primi_vicini(double *r, double *r_prim_vic, double distanza_interaz, int N_mol, int j){
    for(int i=0; i<3*N_mol; i+=3){
        r_prim_vic[i]=r[i]+distanza_interaz*round(abs(r[j]-r[i])/distanza_interaz);
        r_prim_vic[i+1]=r[i+1]+distanza_interaz*round(abs(r[j+1]-r[i+1])/distanza_interaz);
        r_prim_vic[i+2]=r[i+2]+distanza_interaz*round(abs(r[j+2]-r[i+2])/distanza_interaz);
    }
}

#endif /* funzioni_utili_h */
