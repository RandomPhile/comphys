#include <iostream>
#include <cmath>

using namespace std;

/*Due metodi per calcolare la derivata:

f'(x0,h) = [f(x0+h) - f(x0)] / h       - f''(chi)*h/2
f'(x0,h) = [f(x0+h) - f(x0-h)] / (2h)  - f'''(chi)*h^2/6

*/ 

void lista_n(double x_min, double x_max, int n, double *x) {
	double x_step = (x_max-x_min)*(1.0/n);
	for (int j = 0; j <= n; ++j){x[j] = x_min + x_step*j;}
	return;
}

int main() {
	double x0 = 1.0;
	double h_min = 0.01, h_max = 1;
	int n = 100;

	double h[n],err1[n],err2[n];

	lista_n(h_min,h_max,n,h);
	
	for (int i = 0; i <= n; ++i) {
		err1[i] = fabs(exp(x0) - ((exp(x0 + h[i]) - exp(x0))/h[i]));
		err2[i] = fabs(exp(x0) - ((exp(x0 + h[i]) - exp(x0 - h[i]))/(2.0*h[i])));
		cout << h[i] << "\t" << err1[i] << "\t" << err2[i] << endl;
	}
}
