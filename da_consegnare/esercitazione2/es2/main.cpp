#include <iostream>
#include <cmath>
#include <fstream>
#include "funzioni.h"
#define _USE_MATH_DEFINES
using namespace std;
ofstream dati;

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
    printf("v = %f\n\n", v);

    dati << "Metodo bisezione\n" << endl;
    e = bisezione(0, 2, delta);
    printf("e = %f\tE = %f\tSteps = %d\n", e, -e * lam, contatore);

    dati << "\n\nMetodo secante\n" << endl;
    contatore = 0;
    e = secante(1.9, 2, delta);
    printf("e = %f\tE = %f\tSteps = %d\n", e, -e * lam, contatore);

    dati << "\n\nMetodo Newton-Raphson\n" << endl;
    contatore = 0;
    e = NR(2, 4);
    printf("e = %f\tE = %f\tSteps = %d\n", e, -e * lam, contatore);

    /*
    Calcolo raggio quadratico medio usando Simpson
    */
    double E = -e * lam;
    double k = sqrt(2 * muc2 * (V0 + E) / pow(hc, 2));
    double q = sqrt(2 * muc2 * fabs(E) / pow(hc, 2));
    printf("\n\nq=%.7f\tk=%.7f\tR=%.7f\n",q,k,R);





    // printf("\nkR = %f\n\n", k*R);
    // for (int i = 0; i < 2; ++i) {
    //     printf("%f <= kR <= %f\n", (0.5+i)*M_PI, (1+i)*M_PI);
    // }

    double args[] = {0, R, k, 38};
    if ((int)args[3] % 2 != 0) {args[3]++;}
    
    double I1 = simpson(f2, args);
    double I3 = simpson(f1, args);
    
    args[2] = 2 * q;
    
    double I2 = exp(2 * q * R) * pow(sin(k * R), 2) * (2 / pow(2 * q, 3) - simpson(f4, args));
    double I4 = exp(2 * q * R) * pow(sin(k * R), 2) * (1 / (2 * q) - simpson(f3, args));

    printf("I1=%f\tI2=%f\tI3=%f\tI4=%f\n",I1,I2,I3,I4);
    double r2 = (I1 + I2) / (I3 + I4);
    printf("\n<r^2> = %f\n", r2);
    dati.close();
    return 0;
}


