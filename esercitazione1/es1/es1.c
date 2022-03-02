#include <stdio.h>
#include <math.h>

int main(){
	for (float x = 1.f; x > 0; x=x/2){

		if (1.f + x == 1.f){
			printf("Epsilon float  = %.4e\n",x);
			break;
		}
	}
	for (double x = 1.; x > 0; x=x/2){

		if (1. + x == 1.){
			printf("Epsilon double = %.4e\n",x);
			break;
		}
	}
	for (long double x = 1.l; x > 0; x=x/2){

		if (1. + x == 1.){
			printf("Epsilon long   = %.4Le\n",x);
			break;
		}
	}
	return 0;
}