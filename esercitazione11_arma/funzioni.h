#ifndef funzioni_h
#define funzioni_h


#include "header.h"

struct vec {
    double x; double y; double z;

    double mod() {
        return sqrt(x * x + y * y + z * z);
    }
    void uguale(double val) {
        x = val; y = val; z = val;
    }
    void piu(double val) {
        x += val; y += val; z += val;
    }
    void per(double val) {
        x *= val; y *= val; z *= val;
    }
    void print(){
        cout << "x=" <<x << ", y=" << y << ", z=" << z << endl;
    }
};

double pow1(double base, int esp) {//funzione esponenziale creata per non usare pow
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}
double min(double a, double b){
    if(a<=b){
        return a;
    }
    else{
        return b;
    }
}
double Fabs(int x){
    if(x<0){
        return -x;
    }
    else{
        return x;
    }
}

#endif 
