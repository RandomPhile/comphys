#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
double log1_(double x_);//definisco due funzioni per avere la distinzione tra float e double per la funzione con la correzione per la stabilità
float log1(float x);

int main() {
    ofstream dati1, dati2;
    dati1.open("dati1.dat"); dati2.open("dati2.dat");
    double x1_ = 1.0, x2_ = 0.1, err_rel1, err_rel2, err_rel3, err_rel4;
    float x1 = (float) x1_, x2 = (float) x2_;

    for (int i = 0; i < 40; i++) {
        err_rel1 = fabs((log1p(x1_) - log1pf(x1)) / log1p(x1_));
        err_rel2 = fabs((log1p(x2_) - log1pf(x2)) / log1p(x2_));

        err_rel3 = fabs((log1_(x1_) - log1(x1)) / log1_(x1_));
        err_rel4 = fabs((log1_(x2_) - log1(x2)) / log1_(x2_));
        
        dati1 << x1_ << "\t" << err_rel1 << "\t" << err_rel3 << endl;
        dati2 << x2_ << "\t" << err_rel2 << "\t" << err_rel4 << endl;
        
        x1 /= 2.f; x2 /= 2.f;
        x1_ /= 2.0; x2_ /= 2.0;
    }
    dati1.close(); dati2.close(); return 0;
}


double log1_(double x_) {
    if (x_ + 1.0 == 1.0) {
        return x_;
    } else {
        double y_ = (log(1.0 + x_) / ((1.0 + x_) - 1.0));
        return x_ * y_;
    }
}

float log1(float x) { //funzioni definite come secondo l'esercizio della scheda
    if (x + 1.f  == 1.f) {
        return x;
    } else {
        float y = (logf(1.f + x) / ((1.f + x) - 1.f));
        return x * y;
    }
}