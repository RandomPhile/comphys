#include <iostream>

using namespace std;

/*Due metodi per calcolare la derivata:
Differenza in avanti
f'(x0,h) = [f(x0+h) - f(x0)] / h       - f''(chi)*h/2

Differenza centrale
f'(x0,h) = [f(x0+h) - f(x0-h)] / (2h)  - f'''(chi)*h^2/6

*/ 

//questa funzione crea un array di n elementi da un valore minimo ad uno max
//in realtà è un void (non ritorna nulla) perchè modifica un array x dato come argomento
//(più avanti potremo sostituirla con una libreria tipo <vectors>)
void lista_n_double(double x_min, double x_max, int n, double *x) {
	double x_step = (x_max-x_min)*(1.0/n);
	for (int i = 0; i < n; ++i){x[i] = x_min + x_step*i;}
	return;
}

int main() {
	//*** parametri da cambiare eventualmente:
	int n = 100;
	double x0 = 1.0;
	double h[n], h_min = 0.01, h_max = 1;
	//***

	lista_n_double(h_min,h_max,n,h);//assegno valori all'array h[]
	double err1[n], err2[n];
	
	for (int i = 0; i < n; ++i) {
		err1[i] = fabs(exp(x0) - ((exp(x0 + h[i]) - exp(x0))/h[i]));//primo metodo di calcolo derivata
		err2[i] = fabs(exp(x0) - ((exp(x0 + h[i]) - exp(x0 - h[i]))/(2.0*h[i])));//secondo metodo
		cout << h[i] << "\t" << err1[i] << "\t" << err2[i] << endl;//stampa 3 colonne: h err1 err2
	}
}
