/*
L'idea è quella di stimare l'errore macchina "epsilon" cercando il 
valore più grande di x che rispetta l'uguaglianza (1 + x == 1), 
dove x ed 1 possono essere float, double o long double.
*/

#include <iostream>
using namespace std;

int main() {
	for (float x = 1.f; x > 0; x = x / 2) {
		if (1.f + x == 1.f) {
			cout << "Epsilon float  = " << x << endl;
			break;
		}
	}

	for (double x = 1.0; x > 0; x = x / 2) {
		if (1.0 + x == 1.0) {
			cout << "Epsilon double = " << x << endl;
			break;
		}
	}

	for (long double x = 1.l; x > 0; x = x / 2) {
		if (1.l + x == 1.l) {
			cout << "Epsilon long   = " << x << endl;
			break;
		}
	}

	return 0;
}
