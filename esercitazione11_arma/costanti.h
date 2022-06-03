#include <iostream>
#include <fstream>
#include <cmath>

#define ARMA_NO_DEBUG//non fa controlli e va piu veloce, da togliere solo se il codice funziona gia perfettamente

#include <armadillo>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace arma;



int caso_min = 9;//mettere -1 per avere P(rho)




/*** variabili globali ***/
//CC, BCC, FCC
int M = 1; //1,2,4
int N = M * pow(6, 3); //numero di particelle
int N_t;//numero passi simulazione, vedi rigo 61 per il valore in base al tempo di equilibrazione stimato

int N_b=20;//numero bin
double T_req = 1.1;//temperatura adimensionale

string coord_path = "out/coord.dat";
string dati_path = "out/dati.dat";
string gnuplot_path = "out/gnuplot.dat";
string oss_path = "out/oss.dat";

int numero_proposti=0;
int numero_accettati=0;

#define LOG(x) cout<<x<<endl;