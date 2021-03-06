#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

int main() {
    ofstream dati1, dati2;
    dati1.open("dati1.dat"); dati2.open("dati2.dat");
    //*** parametri da cambiare eventualmente:
    const int n = 1000;
    double a_ = 1.0, c_ = 1.0;//scelti arbitrariamente, con b^2>>4ac
    double b_min_ = 2, b_max_ = 5;
    //***
    double x1_, x2_, x1s_, x2s_, err_rel1, err_rel2, err_rel1s, err_rel2s;
    float x1, x2, x1s, x2s;
    float a = (float) a_, b, c = (float) c_;
    
    for (double b_ = pow(10,b_min_); fabs(b_) < pow(10,b_max_); b_ *= pow(10,(b_max_-b_min_)/n)) {
        //per valori di b<0:
        //b_ = -fabs(b_);
        b = (float) b_;

        //check perchè ho una radice quadrata
        if ((powf(b,2) - 4.f*a*c < 0) || (pow(b_,2) - 4.0*a_*c_ < 0)) {return 1;}
        
        x1 = (-b + sqrt(powf(b,2) - 4.f*a*c)) / (2.f*a);
        x2 = (-b - sqrt(powf(b,2) - 4.f*a*c)) / (2.f*a);
        x1_ = (-b_ + sqrt(pow(b_,2) - 4.0*a_*c_)) / (2.0*a_);
        x2_ = (-b_ - sqrt(pow(b_,2) - 4.0*a_*c_)) / (2.0*a_);
        
        err_rel1 = fabs((x1_ - x1) / x1_);
        err_rel2 = fabs((x2_ - x2) / x2_);
        
        //soluzione alternativa:
        x1s = (2.f*c) / (-b - sqrt(powf(b,2) - 4.f*a*c));
        x2s = (2.f*c) / (-b + sqrt(powf(b,2) - 4.f*a*c));
        
        x1s_ = (2.0*c_) / (-b_ - sqrt(pow(b_,2) - 4.0*a_*c_));
        x2s_ = (2.0*c_) / (-b_ + sqrt(pow(b_,2) - 4.0*a_*c_));
        
        err_rel1s = fabs((x1s_ - x1s) / x1s_);
        err_rel2s = fabs((x2s_ - x2s) / x2s_);

        dati1 << fabs(b_) << "\t" << err_rel1 << "\t" << err_rel2 << endl;
        dati2 << fabs(b_) << "\t" << err_rel1s << "\t" << err_rel2s << endl;
    }
    dati1.close(); dati2.close();
    return 0;
}
