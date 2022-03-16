#include <iostream>

//Per evitare di scrivere std:: ogni volta prima di cout, endl, ecc.. usiamo:
using namespace std;

int main() {

	for (float x = 1.f; x > 0; x=x/2){
		if (1.f + x == 1.f){//uso == per vedere quando trovo un numero che sommato a 1 non cambi il risultato, quindi l'epsilon di macchina
			cout<<"Epsilon float  = "<<x<<endl;
			break;
		}
	}

	for (double x = 1.0; x > 0; x=x/2){
		if (1.0 + x == 1.0){
			cout<<"Epsilon double = "<<x<<endl;
			break;
		}
	}
	
	for (long double x = 1.l; x > 0; x=x/2){
		if (1.l + x == 1.l){
			cout<<"Epsilon long   = "<<x<<endl;
			break;
		}
	}
	
    return 0;
}
