#include <stdio.h>
#include <math.h>

int main(){
	float x0 = 1.f;
	int n = 100;
	double h[n];
	double err[n];
//ciao
	for (int j = 0; j < n; ++j){
		h[j] = (j+1)*(0.000001l/n);
	}
	

	for (int i = 0; i < n; ++i) {
		err[i] = fabs(exp(x0) - ((exp(x0 + h[i]) - exp(x0))/h[i]));
		printf("%10.8f\t%10.8f\n",h[i],err[i]);
	}
}
