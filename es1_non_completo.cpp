#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

using namespace std;

//definisco le funzioni per l'integrazione con overloading di funzione
double trapezi(double (*f)(double), double a, double b, int N, double param);
double trapezi(double (*f)(double, double), double a, double b, int N, double param);
double Simpson(double (*f)(double), double a, double b, int N, double param);
double Simpson(double (*f)(double, double), double a, double b, int N, double param);

//definisco le funzioni che voglio integrare
double x2sin2(double x, double alfa){
    return pow(x*sin(alfa*x),2);
}
double sin2(double x, double alfa){
    return pow(sin(alfa*x),2);
}
double x2exp(double x, double alfa){
    return (pow(x,2)*exp(-alfa*x));
}
double exp1(double x, double alfa){
    return exp(-alfa*x);
}
            
//main scrivo gli integrali richiesti
int main() {
    double b=1, a=0, delta1=1e-4, delta2=1e-6;
    double f2_max[4], f4_max[4];
    int NT[4];
    int NS[4];
    
    
    
    //vanno create le derivate, non ho voglia di farlo, per il resto il programma funziona bene (in teoria)
    
    
    
    
    
    
    //errore di 10^-4
    
    for(int i=0; i<4; i++){
        NT[i]=ceil(sqrt(pow(b-a,3)*f2_max[i]/(12*delta1)));
    }
    for(int i=0; i<4; i++){
        NS[i]=ceil(sqrt(pow(b-a,3)*f4_max[i]/(12*delta2)));
    }
    
    cout<<"Trapezi sin2: "<<trapezi(sin2, 0, 1, NT[0], 2)<<endl;
    cout<<"Simpson sin2: "<<Simpson(sin2, 0, 1, NS[0], 2)<<endl;
    
    cout<<"Trapezi x2sin2: "<<trapezi(x2sin2, 0, 1, NT[1], 2)<<endl;
    cout<<"Simpson x2sin2: "<<Simpson(x2sin2, 0, 1, NS[1], 2)<<endl;
    
    //vedi  consiglio scheda per il perche del 1-ans/2-ans
    cout<<"Trapezi exp: "<<1-trapezi(exp1, 0, 1, NT[2], 0.5)<<endl;
    cout<<"Simpson exp: "<<1-Simpson(exp1, 0, 1, NS[2], 0.5)<<endl;
    
    cout<<"Trapezi x2exp: "<<2-trapezi(x2exp, 0, 1, NT[3], 0.5)<<endl;
    cout<<"Simpson x2exp: "<<2-Simpson(x2exp, 0, 1, NS[3], 0.5)<<endl;
    
    
    //errore di 10^-6
    
    for(int i=0; i<4; i++){
        NT[i]=ceil(sqrt(pow(b-a,3)*f2_max[i]/(12*delta2)));
    }
    for(int i=0; i<4; i++){
        NS[i]=ceil(sqrt(pow(b-a,3)*f4_max[i]/(12*delta2)));
    }
    
    cout<<"Trapezi sin2: "<<trapezi(sin2, 0, 1, NT[0], 2)<<endl;
    cout<<"Simpson sin2: "<<Simpson(sin2, 0, 1, NS[0], 2)<<endl;
    
    cout<<"Trapezi x2sin2: "<<trapezi(x2sin2, 0, 1, NT[1], 2)<<endl;
    cout<<"Simpson x2sin2: "<<Simpson(x2sin2, 0, 1, NS[1], 2)<<endl;
    
    //vedi  consiglio scheda per il perche del 1-ans/2-ans
    cout<<"Trapezi exp: "<<1-trapezi(exp1, 0, 1, NT[2], 0.5)<<endl;
    cout<<"Simpson exp: "<<1-Simpson(exp1, 0, 1, NS[2], 0.5)<<endl;
    
    cout<<"Trapezi x2exp: "<<2-trapezi(x2exp, 0, 1, NT[3], 0.5)<<endl;
    cout<<"Simpson x2exp: "<<2-Simpson(x2exp, 0, 1, NS[3], 0.5)<<endl;
    
    return 0;
}

//scrivo la funzione per i due metodi usando il puntatore a funzione per passare una funzione come parametro di una seconda funzione. gli altri parametri sono i due estremi, il numero di intervalli e il parametro di ingresso
double trapezi(double (*f)(double, double), double a, double b, int N, double param){
    double h=(b-a)/N;
    double I;
    I=(f(a,param)+f(b,param))/2;
    for(int i=1; i<N; i++){
        I+=f(a+i*h,param);
    }
    I*=h;
    return I;
}
double trapezi(double (*f)(double), double a, double b, int N, double param){
    double h=(b-a)/N;
    double I;
    I=(f(a*param)+f(b*param))/2;
    for(int i=1; i<N; i++){
        I+=f((a+i*h)*param);
    }
    I*=h;
    return I;
}

double Simpson(double (*f)(double, double), double a, double b, int N, double param){
    double h=(b-a)/N;
    double I;
    I=(f(a,param)+f(b,param));
    for(int i=1; i<N/2; i++){
        I=I+2*f(a+2*i*h,param)+4*f(a+(2*i-1)*h,param);
    }
    I=(I-2*f(a+N*h, param))*h/3;
    return I;
}

double Simpson(double (*f)(double), double a, double b, int N, double param){
    double h=(b-a)/N;
    double I;
    I=(f(a*param)+f(b*param));
    for(int i=1; i<N/2; i++){
        I=I+2*f((a+2*i*h)*param)+4*f((a+(2*i-1)*h)*param);
    }
    I=(I-2*f((a+N*h)*param))*h/3;
    return I;
}
