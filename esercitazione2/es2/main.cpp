#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
ofstream dati;

double f(double x);
double bisezione(double a, double b, double d);
double secante(double a, double b, double d);
double NR(double a, int max);//Newton-Raphson

double v;
int contatore = 0;

int main() {
    dati.open("dati.dat");
    double e;
    double delta = 1e-6;
    double hc = 197.327;//MeVfm
    double R = 1.93;//fm
    double V0 = 38.5;//MeV
    //V0 = 27.4;
    double Mnc2 = 939.565;//MeV
    double Mpc2 = 938.272;//MeV

    double muc2 = Mnc2 * Mpc2 / (Mnc2 + Mpc2);
    double lam = pow(hc, 2) / (2 * muc2 * pow(R, 2));
    v = 2 * V0 * muc2 * pow(R, 2) / pow(hc, 2);
    printf("v = %f\n",v);
    


    dati << "Metodo bisezione\n" << endl;
    e = bisezione(0, 2, delta);
    printf("e = %f\tE = %f\tSteps = %d\n", e, -e * lam, contatore);

    dati << "\n\nMetodo secante\n" << endl;
    contatore = 0;
    e = secante(1.9, 2, delta);
    printf("e = %f\tE = %f\tSteps = %d\n", e, -e * lam, contatore);
    
    dati << "\n\nMetodo Newton-Raphson\n" << endl;
    contatore = 0;
    e = NR(2,20);
    printf("e = %f\tE = %f\tSteps = %d\n", e, -e * lam, contatore);
    
    dati.close();
}

double f(double x) {
    return 1 / tan(sqrt(v - x)) + sqrt(x / (v - x));
}

double f_prime(double x) {
    return 1 / (2 * pow(sin(sqrt(v - x)), 2) * sqrt(v - x)) + v / (2 * sqrt(x / (v - x)) * pow((x - v), 2));
}

double bisezione(double a, double b, double d) {
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
        c = b - f(b) * (b - a) / (f(b) - f(a));
        a = b;
        b = c;
        dati << contatore << '\t' << c << endl;
        contatore++;
    }
    return c;
}

double NR(double a, int max) {
    while (contatore <= max) {
        a = a - f(a)/f_prime(a);
        dati << contatore << '\t' << a << endl;
        contatore++;
    }
    return a;
}
