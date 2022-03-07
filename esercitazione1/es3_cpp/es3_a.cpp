#include <iostream>

using namespace std;

//solite funzioni per creare array in un range (vedi es2_a)
//questa volta anche una coi float
void lista_n_float(float x_min, float x_max, int n, float *x) {
    float x_step = (x_max-x_min)*(1.f/n);
    for (int i = 0; i < n; ++i){x[i] = x_min + x_step*i;}
    return;
}
void lista_n_double(double x_min, double x_max, int n, double *x) {
    double x_step = (x_max-x_min)*(1.0/n);
    for (int i = 0; i < n; ++i){x[i] = x_min + x_step*i;}
    return;
}

/*
scelgo i valori di "c" tale che il discriminante sia sempre positivo:
b^2 - 4ac >= 0 
c <= b^2/4a
*/

int main() {
    //*** parametri da cambiare eventualmente:
    int n = 100;
    float  a  = 1e-6, b_min  = 1e4, b_max  = 1e5;
    double a_ = 1e-6, b_min_ = 1e4, b_max_ = 1e5;
    //***

    //float
    float c, b[n], x1[n], x2[n];
    
    if (b_min >= 0) {//modo brutto brutto per assicurarmi che discriminante sempre >= 0 e ho 2 soluzioni reali
        c = (b_min*b_min/(4.f*a)) - 0.1f;
    } else {
        c = -0.1f;
    }
    
    lista_n_float(b_min,b_max,n,b);
    
    //double (probabilmente sono più double di quanti servono, comunque le variabili che finiscono con trattino basso _ sono double)
    double c_, b_[n], x1_[n], x2_[n];

    if (b_min >= 0) {
        c_ = (b_min_*b_min_/(4.0*a_)) - 0.1;
    } else {
        c_ = -0.1;
    }
    
    lista_n_double(b_min_,b_max_,n,b_);

    //questo serve per il programma di gnuplot (vedi run_es3_a.sh)
    cout << "# Soluzione x1 (index 0)" << endl;

    for (int i = 0; i < n; ++i) {
        //float
        if (b[i]*b[i] - 4.f*a*c < 10e-8) {//check perchè ho una radice quadrata
            printf("Discriminante negativo = %f\n", b[i]*b[i] - 4.f*a*c);//non dovrebbe mai succedere
        } else {
            x1[i] = (-b[i] + sqrt(b[i]*b[i] - 4.f*a*c))/(2.f*a);
            x2[i] = (-b[i] - sqrt(b[i]*b[i] - 4.f*a*c))/(2.f*a);
        }

        //double
        if (b_[i]*b_[i] - 4.f*a_*c_ < 10e-16) {//check perchè ho una radice quadrata
            printf("Discriminante negativo = %f\n", b_[i]*b_[i] - 4.f*a_*c_);//non dovrebbe mai succedere
        } else {
            x1_[i] = (-b_[i] + sqrt(b_[i]*b_[i] - 4.0*a_*c_))/(2.0*a_);
            x2_[i] = (-b_[i] - sqrt(b_[i]*b_[i] - 4.0*a_*c_))/(2.0*a_);
        }

        //stampo in 3 colonne b sol1_float sol1_double
        cout << b[i] << "\t" << x1[i] << "\t" << x1_[i] << endl;
    }
    
    cout << "\n\n# Soluzione x2 (index 1)" << endl;

    for (int i = 0; i <  n; ++i) {
        //stampo in 3 colonne b sol2_float sol2_double
        cout << b[i] << "\t" << x2[i] << "\t" << x2_[i] << endl;
    }
    return 0;
}
