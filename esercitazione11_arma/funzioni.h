#ifndef funzioni_h
#define funzioni_h


#include "header.h"

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
#endif 
