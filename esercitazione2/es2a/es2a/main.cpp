#include <iostream>
#include <cmath>
#include <fstream>
#define _USE_MATH_DEFINES
using namespace std;
ofstream dati;

double f1(double x);
double fA(double A, double r);
double bisezione(double (*f)(double), double a, double b, double d);
double bisezione(double (*f)(double, double), double a, double b, double d, double A);
double secante(double a, double b, double d);
double NR(double a, int max);//Newton-Raphson
double psi2(double r, double A);
double r2psi2(double r, double A);
double max_of_function(double (*f)(double, double, int), int der, double a, double b, double param, int N_der);
double simpson(double (*f)(double, double), double a, double b, int N, double param);


double v;
double R;
double k;
double q;
int contatore = 0;

int main() {
    dati.open("dati.dat");
    double e;
    double delta = 1e-6;
    double hc = 197.327;//MeVfm
    R = 1.93;//fm
    double V0 = 38.5;//MeV
    //V0 = 27.4;
    double Mnc2 = 939.565;//MeV
    double Mpc2 = 938.272;//MeV

    double muc2 = Mnc2 * Mpc2 / (Mnc2 + Mpc2);
    double lam = pow(hc, 2) / (2 * muc2 * pow(R, 2));
    v = 2 * V0 * muc2 * pow(R, 2) / pow(hc, 2);
    printf("v = %f\n",v);
    
    
    dati << "Metodo bisezione\n" << endl;
    e = bisezione(f1,0, 2, delta);
    printf("e = %f\tE = %f\tSteps Bis = %d\n", e, -e * lam, contatore);

    dati << "\n\nMetodo secante\n" << endl;
    contatore = 0;
    e = secante(1.9, 2, delta);
    printf("e = %f\tE = %f\tSteps Sec = %d\n", e, -e * lam, contatore);
    
    dati << "\n\nMetodo Newton-Raphson\n" << endl;
    contatore = 0;
    e = NR(2,20);
    printf("e = %f\tE = %f\tSteps NR = %d\n", e, -e * lam, contatore);
    
    cout<<"\n\n";
    

    //devo ora trovare il raggio quadratico medio usando simpson. copio il codice dall'esercizio uno e definisco la psi da usare per l'integrale
    
    double r_2M;
    //pongo k,R,q variabili esterne per poter usare il programma trovato in precedenza
    k=sqrt(2*muc2*(V0-e*lam)/pow(hc,2));
    q=sqrt(2*muc2*e*lam/pow(hc,2));
    printf("k: %f\tq: %f\n",k,q);
    //trovo il valore di A per la quale io abbia psi normalizzata con la bisezione
    double A_min=2, A_max=4, A=0;
    while (abs(2 * (A_max - A_min) / (A_max + A_min)) > delta) {
        A = (A_max + A_min) / 2;
        double I1 = simpson(psi2, 1e-3, 1e6, 1e7, A_max);
        double I2 = simpson(psi2, 1e-3, 1e6, 1e7, A);
        
        if ((I1-1)*(I2-1) < 0){
            A_min = A;
           
        }
        else {
            A_max = A;
        }
        
    }
    cout<<A<<endl;
    r_2M=simpson(r2psi2, 1e-3, 1e5, 1e7, A);
    cout<<r_2M<<endl;
    dati.close();
}

double f1(double x) {
    return 1 / tan(sqrt(v - x)) + sqrt(x / (v - x));
}
double fA(double A, double r){
    return A*(sin(k*r)+sin(k*R)*exp(q*(R-r)));
}

double f_prime(double x) {
    return 1 / (2 * pow(sin(sqrt(v - x)), 2) * sqrt(v - x)) + v / (2 * sqrt(x / (v - x)) * pow((x - v), 2));
}

double bisezione(double (*f)(double), double a, double b, double d){
    double c;
    while (abs(2 * (b - a) / (a + b)) > d) {
        c = (a + b) / 2;
        if (f(a)*f(c) < 0) {
            b = c;
        } else {
            a = c;
        }
        dati << contatore << '\t' << c << endl;
        contatore++;
    }
    return c;
}

double secante(double a, double b, double d) {
    double c;
    while (abs((b - a) / b) > d) {
        c = b - f1(b) * (b - a) / (f1(b) - f1(a));
        a = b;
        b = c;
        dati << contatore << '\t' << c << endl;
        contatore++;
    }
    return c;
}

double NR(double a, int max) {
    while (contatore <= max) {
        a = a - f1(a)/f_prime(a);
        dati << contatore << '\t' << a << endl;
        contatore++;
    }
    return a;
}
double simpson(double (*f)(double, double), double a, double b, int N, double param) {
    double h = (b - a) / N;
    //cout<<h<<endl;
    double I;
    I = (f(a, param) + f(b, param));
    for (int n = 1; n <= N / 2; n++) {
        I += 2 * f(a + 2 * n * h, param) + 4 * f(a + (2 * n - 1) * h, param);
    }
    //tolgo l'elemento N/2 della prima sommatoria
    I = (I - 2 * f(b, param)) * h / 3;
    return I;
}
double max_of_function(double (*f)(double, double, int), int der, double a, double b, double param, int N_der) {
    float temp[N_der];
    double h_der;
    for (int n = 0; n < N_der; n++) {
        h_der = a + n * b / (1.0 * N_der);
        temp[n] = abs(f(h_der, param, der));
    }
    return *max_element(temp, temp + N_der);
}

double psi2(double r, double A){
    if(r<=R){
        double y=A*sin(k*r)/(r*sqrt(4*M_PI));
        return y*y;
    }
    else{
        double y=A*sin(k*R)*exp(q*(R-r))/(r*sqrt(4*M_PI));
        return y*y;
    }
}
double r2psi2(double r, double A){
    if(r<=R){
        double y=r*A*sin(k*r)/(r*sqrt(4*M_PI));
        return y*y;
    }
    else{
        double y=r*A*sin(k*R)*exp(q*(R-r))/(r*sqrt(4*M_PI));
        return y*y;
    }
}
