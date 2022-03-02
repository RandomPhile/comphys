#include <stdio.h>
#include <math.h>

int main(){
	float x0 = 1;
	int n = 1000;
	float h[n];
	float err[n];

	for (int j = 0; j < n; ++j){
		h[j] = (j+1)*(0.01/n);
		//printf("%d\t%d\t%f\n",j,n,h[j]);
	}
	

	for (int i = 0; i < n; ++i) {
		err[i] = fabs(exp(x0) - ((exp(x0 + h[i]) - exp(x0))/h[i]));
		printf("%10.8f\t%10.8f\n",h[i],err[i]);
	}
}