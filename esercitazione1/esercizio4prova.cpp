#include <iostream>
#include <cmath>

double log1(double x);//definisco due funzioni per avere la distinzione tra float e double per la funzione con la correzione per la stabilità
float log21(float x);

using namespace std;

int main() {
    float x0=1.f;
    double x1=1.l, err_rel, err_rel1;
    
    for(int i=0; i<50; i++){
        err_rel=fabs((log(1.l+x1)-log(1.f+x0))/log(1.l+x1));//osservo gli errori relativi per i vari valori di 1+x del caso della funzione log del compilatore e per la funzione creata stabile
        err_rel1=fabs((log1(1.l+x1)-log21(1.f+x0))/log1(1.l+x1));
        cout<<x0<<"\t"<<err_rel<<"\t"<<err_rel1<<endl;//stampo a video il valore di x0 e gli errori relativi per le due funzioni
        x0=x0/2.f;
        x1=x1/2.l;
    }
    
    return 0;
}

double log1(double x){//funzioni definite come secondo l'esercizio della scheda
    if(1.l+x==1.l){
        return x;
    }
    else{
        double y;
        y=x*(log(1.l+x)/((1.l+x)-1.l));
        return y;
    }
}
float log21(float x){
    if(1.f+x==1.f){
        return x;
    }
    else{
        double y;
        y=x*(log(1.f+x)/((1.f+x)-1.f));
        return y;
    }
}
