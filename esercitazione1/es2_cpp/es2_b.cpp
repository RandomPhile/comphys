#include <iostream>
#include <cmath>

using namespace std;

/*Errore nel fare l'operazione @:

(a @ b) -> (a @ b)*(1+E/2)

ed E può essere anche negativo
*/ 

//stessa funzione per gli array (vedi es2_a)
void lista_n_double(double x_min, double x_max, int n, double *x) {
	double x_step = (x_max-x_min)*(1.0/n);
	for (int i = 0; i < n; ++i){x[i] = x_min + x_step*i;}
	return;
}

int main(){
	//*** parametri da cambiare eventualmente:
	int n = 100;
	double x0 = 1.0, E = 0.1;
	double h[n], h_min = 0.01, h_max = 2;
	//***
	
	double E_ = 1.0 + E/2.0;//per comodità
	lista_n_double(h_min,h_max,n,h);//assegno array
	double err1[n], err2[n], err3[n];
	
	for (int i = 0; i < n; ++i) {
		/*prima formula per la derivata con l'errore su 3 operazioni:
		# x0+h
		# f(x0+h) - f(x0)
		# [f(x0+h) - f(x0)] / h
		*/
		err1[i] = fabs(exp(x0) - ((exp((x0 + h[i])*E_) - exp(x0))*(E_*E_)/h[i]));

		/*seconda formula per la derivata con l'errore su 4 operazioni:
		# x0+h
		# x0-h
		# f(x0+h) - f(x0-h)
		# [f(x0+h) - f(x0-h)] / 2h
		*/
		err2[i] = fabs(exp(x0) - ((exp((x0 + h[i])*E_) - exp((x0 - h[i])*E_))*(E_*E_)/(2.0*h[i])));
		
		/*seconda formula per la derivata con l'errore su 5 operazioni (così facendo si semplifica un E_:
		# x0+h
		# x0-h
		# f(x0+h) - f(x0-h)
		# 2.0*h
		# [f(x0+h) - f(x0-h)] / 2h
		*/
		err3[i] = fabs(exp(x0) - ((exp((x0 + h[i])*E_) - exp((x0 - h[i])*E_))*E_/(2.0*h[i])));
		
		cout << h[i] << "\t" << err1[i] << "\t" << err2[i] << "\t" << err3[i] << endl;
	}
}
