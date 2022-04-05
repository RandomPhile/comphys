/*
Usiamo due metodi per calcolare la derivata:
Differenza in avanti
f'(x0,h) = [f(x0+h) - f(x0)] / h       - f''(chi)*h/2

Differenza centrale
f'(x0,h) = [f(x0+h) - f(x0-h)] / (2h)  - f'''(chi)*h^2/6

Per stimare l'errore in modulo facciamo la differenza
con la derivata effettiva e^x
Ripetiamo l'operazione con variabili float e double
*/

#include <iostream>
#include <cmath>
#include <fstream>//per creare e scrivere sul file dati.dat
using namespace std;
ofstream dati;

int main() {
	dati.open("dati.dat");
	//*** parametri da cambiare eventualmente:
	int n = 120;
	double x0 = 1.0;
	double he_min = -18, he_max = 0;
	//***

	double err1, err2;
	int i = 0;
	for (double h = pow(10, he_min); h < pow(10, he_max); h *= pow(10, (he_max - he_min) / n)) {
		err1 = fabs(exp(x0) - ((exp(x0 + h) - exp(x0)) / h)); //primo metodo di calcolo derivata
		err2 = fabs(exp(x0) - ((exp(x0 + h) - exp(x0 - h)) / (2.0 * h))); //secondo metodo
		dati << h << "\t" << err1 << "\t" << err2 << endl; //scrivo 3 colonne: h err1 err2
		i++;
	}
	dati << "\n\n" << endl; //nuovo blocco per gnuplot

	float x0f = (float) x0, he_minf = (float) he_min, he_maxf = (float) he_max, err1f, err2f;
	i = 0;
	for (float h = powf(10, he_minf); h < powf(10, he_maxf); h *= powf(10, (he_maxf - he_minf) / n)) {
		err1f = fabs(exp(x0f) - ((exp(x0f + h) - exp(x0f)) / h)); //primo metodo di calcolo derivata
		err2f = fabs(exp(x0f) - ((exp(x0f + h) - exp(x0f - h)) / (2.0 * h))); //secondo metodo
		dati << h << "\t" << err1f << "\t" << err2f << endl; //scrivo 3 colonne: h err1f err2f
		i++;
	}
	
	dati.close();
}