#include <iostream>
#include <cmath>

using namespace std;

/*Errore nel fare l'operazione @
(a @ b) -> (a @ b)*(1+E/2)
ed E può essere anche negativo
*/ 

void lista_n_double(double x_min, double x_max, int n, double *x) {
	double x_step = (x_max-x_min)*(1.0/n);
	for (int i = 0; i < n; ++i){x[i] = x_min + x_step*i;}
	return;
}

int main(){
	double x0 = 1.0, E = 0.01;
	double E_ = 1.0 + E/2.0;
	int n = 100;
	double h[n], h_min = 0.01, h_max = 2;
	lista_n_double(h_min,h_max,n,h);
	double err1[n], err2[n], err3[n];
	
	for (int i = 0; i < n; ++i) {
		err1[i] = fabs(exp(x0) - ((exp((x0 + h[i])*E_) - exp(x0))*(E_*E_)/h[i]));
		err2[i] = fabs(exp(x0) - ((exp((x0 + h[i])*E_) - exp((x0 - h[i])*E_))*(E_*E_)/(2.0*h[i])));
		//Se considero anche l'errore in 2*h:
		err3[i] = fabs(exp(x0) - ((exp((x0 + h[i])*E_) - exp((x0 - h[i])*E_))*E_/(2.0*h[i])));
		
		cout << h[i] << "\t" << err1[i] << "\t" << err2[i] << "\t" << err3[i] << endl;
	}
}
