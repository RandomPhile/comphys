#include <iostream>
#include <cmath>
#include <fstream>

#include "funzioni.h"
#include "integrali.h"
#include "zeri.h"
#include "psi.h"

#define _USE_MATH_DEFINES
using namespace std;
ofstream dati;

double f1(double x);//funzione iniziale per trovare E

double v;
double R;
double k;
double q;
int contatore=0;

int main() {
    dati.open("dati.dat");
    
    //definizione variabili utilizzate
    double e;
    double delta = 1e-6;
    double hc = 197.327;//MeVfm
    R = 1.93;//fm
    double V0 = 38.5;//MeV
    //V0 = 27.4;
    double Mnc2 = 939.565;//MeV
    double Mpc2 = 938.272;//MeV

    //calcolo variabili richieste
    double muc2 = Mnc2 * Mpc2 / (Mnc2 + Mpc2);
    double lam = pow(hc, 2) / (2 * muc2 * pow(R, 2));
    v = 2 * V0 * muc2 * pow(R, 2) / pow(hc, 2);
    printf("v = %f\n",v);
    
    //calcolo di E coi vari metodi
    dati << "Metodo bisezione\n" << endl;
    e = bisezione(f1,0, 2, delta);
    printf("e = %f\tE = %f\tSteps Bis = %d\n", e, -e * lam, contatore);

    dati << "\n\nMetodo secante\n" << endl;
    contatore = 0;
    e = secante(f1, 1.9, 2, delta);
    printf("e = %f\tE = %f\tSteps Sec = %d\n", e, -e * lam, contatore);
    
    dati << "\n\nMetodo Newton-Raphson\n" << endl;
    contatore = 0;
    e = NR(f1, 2,20);
    printf("e = %f\tE = %f\tSteps NR = %d\n", e, -e * lam, contatore);
    
    cout<<"\n\n";
    

    //devo ora trovare il raggio quadratico medio usando simpson. copio il codice dall'esercizio uno, lo modifico leggermente e definisco la psi da usare per l'integrale
    
    double r_2M;
    
    //pongo k,R,q variabili esterne per poter usare il programma trovato in precedenza
    k=sqrt(2*muc2*(V0-e*lam)/pow(hc,2));
    q=-k/tan(k*R);
    
 /*   //trovo il valore di A per la quale io abbia psi normalizzata con la bisezione
    double A_min=2, A_max=4, A=0;
    while (fabs(2 * (A_max - A_min) / (A_max + A_min)) > delta) {
        A = (A_max + A_min) / 2;
        double I1 = simpson(psi2, 1e-5, 1e6, 1e7, A_max);
        double I2 = simpson(psi2, 1e-5, 1e6, 1e7, A);
        
        if ((I1-1)*(I2-1) < 0){//uso I1-1 per avere lo zero di funzione trovato con la bisezione
            A_min = A;
           
        }
        else {
            A_max = A;
        }
        
    }*/
    double A=0;
    //cout<<"Il valore di A che normalizza la funzione d'onda è: "<<A<<endl;
    r_2M=simpson(r2psi2, 1e-3, 1e5, 1e7, A)/simpson(psi2, 1e-3, 1e5, 1e7, A);
    cout<<"Il valore del raggio quadratico medio è: "<<r_2M<<" mentre quello esatto (diversa teoria) è: 2.12799"<<endl;
    cout<<"r medio è: "<<simpson(rpsi2, 1e-3, 1e5, 1e7, A)<<endl;
    dati.close();
}




double f1(double x) {
    return 1 / tan(sqrt(v - x)) + sqrt(x / (v - x));
}

