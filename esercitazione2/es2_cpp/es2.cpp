#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
ofstream dati1, dati2;

double f(double x);
double bisezione(double a, double b, double d);

double v;
int contatore = 0;

int main() {
    dati1.open("dati1.dat");
    double delta = 1e-6;
    double hc = 197.327;//MeVfm
    double R = 1.93;//fm
    double V0 = 38.5;//MeV
    double Mnc2 = 939.565;//MeV
    double Mpc2 = 938.272;//MeV

    double muc2 = Mnc2 * Mpc2 / (Mnc2 + Mpc2);
    double lam = pow(hc, 2) / (2 * muc2 * pow(R, 2));
    v = 2 * V0 * muc2 * pow(R, 2) / pow(hc, 2);

    float x_bisezione = bisezione(0, 2, delta);

    double E_bisezione = -x_bisezione * lam;

    printf("e = %f\tE_bisezione = %f\tSteps = %d\n", x_bisezione, E_bisezione, contatore);
    dati1.close();
}

double f(double x) {
    return 1 / tan(sqrt(v - x)) + sqrt(x / (v - x));
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
        dati1 << contatore << '\t' << c << endl;
        contatore++;
    }
    return (a + b) / 2;
}