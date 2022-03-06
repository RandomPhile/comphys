#include <iostream>
#include <cmath>

using namespace std;

/*Due metodi per calcolare la derivata:

f'(x0,h) = [f(x0+h) - f(x0)] / h       - f''(chi)*h/2
f'(x0,h) = [f(x0+h) - f(x0-h)] / (2h)  - f'''(chi)*h^2/6

*/ 

int main(){
	double x0 = 1.0;
	double h_min = 0.01, h_max = 0.1;

	//scegliendo il numero di step:
	int n = 100;
	double h_step = (h_max-h_min)*(1.0/n);
	
	//scegliendo la larghezza degli step:
	//double h_step = 0.001;
	//int n = (h_max-h_min)/h_step;

	double h[n],err1[n],err2[n];

	for (int j = 0; j <= n; ++j){
		h[j] = h_min + h_step*j;
	}
	
	for (int i = 0; i <= n; ++i) {
		err1[i] = fabs(exp(x0) - ((exp(x0 + h[i]) - exp(x0))/h[i]));
		err2[i] = fabs(exp(x0) - ((exp(x0 + h[i]) - exp(x0 - h[i]))/(2.0*h[i])));
		cout << h[i] << "\t" << err1[i] << "\t" << err2[i] << endl;
	}
}
