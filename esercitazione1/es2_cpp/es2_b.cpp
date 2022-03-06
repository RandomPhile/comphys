#include <iostream>
#include <cmath>

using namespace std;

/*Errore nel fare l'operazione @
(a @ b) -> (a @ b)*(1+E/2)
ed E può essere anche negativo
*/ 

int main(){
	double x0 = 1.0, E = 0.01;
	double E_ = 1.0 + E/2.0;
	double h_min = 0.01, h_max = 2;

	//scegliendo il numero di step:
	int n = 10000;
	double h_step = (h_max-h_min)*(1.0/n);
	
	//scegliendo la larghezza degli step:
	//double h_step = 0.001;
	//int n = (h_max-h_min)/h_step;

	double h[n],err1[n],err2[n],err3[n];
	
	for (int j = 0; j <= n; ++j){
		h[j] = h_min + h_step*j;
	}
	
	for (int i = 0; i <= n; ++i) {
		err1[i] = fabs(exp(x0) - ((exp((x0 + h[i])*E_) - exp(x0))*(E_*E_)/h[i]));
		err2[i] = fabs(exp(x0) - ((exp((x0 + h[i])*E_) - exp((x0 - h[i])*E_))*(E_*E_)/(2.0*h[i])));
		//Se considero anche l'errore in 2*h:
		err3[i] = fabs(exp(x0) - ((exp((x0 + h[i])*E_) - exp((x0 - h[i])*E_))*E_/(2.0*h[i])));
		
		cout << h[i] << "\t" << err1[i] << "\t" << err2[i] << "\t" << err3[i] << endl;
	}
}
