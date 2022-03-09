#include <iostream>
#include <cmath>
#include <fstream>

double log1(double x_);//definisco due funzioni per avere la distinzione tra float e double per la funzione con la correzione per la stabilità
float log1_(float x);

using namespace std;

double log1_(double x_) { //funzioni definite come secondo l'esercizio della scheda
    if (x_ + 1.0 == 1.0) {
        return x_;
    } else {
        double y_;
        y_ = (log1p(x_) / (log1p(x_) - 1.0));
        return x_ * y_;
    }
}
double log1(float x) { //funzioni definite come secondo l'esercizio della scheda
    if (x + 1.f  == 1.f) {
        return x;
    } else {
        float y;
        y = (log1pf(x) / (log1pf(x) - 1.f));
        return x * y;
    }
}

int main() {
    ofstream dati1, dati2;
    dati1.open("dati1.dat"); dati2.open("dati2.dat");
    double x1_ = 1.0, x2_ = 0.1, err_rel1, err_rel2, err_rel3, err_rel4;
    float x1 = (float) x1_, x2 = (float) x2_;

    for (int i = 0; i < 100; i++) {
        //cout << log1p(x1_) << "\t" << log1pf(x1) << "\n";
        err_rel1 = fabs((log1p(x1_) - log1pf(x1)) / log1p(x1_));
        err_rel2 = fabs((log1p(x2_) - log1pf(x2)) / log1p(x2_));

        err_rel3 = fabs((log1_(x1_) - log1(x1)) / log1_(x1_));
        err_rel4 = fabs((log1_(x2_) - log1(x2)) / log1_(x2_));
        //err_rel = fabs((log(1.0 + x1) - log(1.f + x0)) / log(1.0 + x1)); //osservo gli errori relativi per i vari valori di 1+x del caso della funzione log del compilatore e per la funzione creata stabile
        //err_rel1 = fabs((log1(1.0 + x1) - log21(1.f + x0)) / log1(1.0 + x1));
        dati1 << x1_ << "\t" << err_rel1 << "\t" << err_rel3 << endl;
        dati2 << x2_ << "\t" << err_rel2 << "\t" << err_rel4 << endl;
        x1 /= 2.f; x2 /= 2.f;
        x1_ /= 2.0; x2_ /= 2.0;
    }
    dati1.close(); dati2.close(); return 0;
}