#ifndef FUNZIONI_H
#define FUNZIONI_H

#include "costanti.h"
#include <iostream>
// #include <fstream>
#include <cmath>
// #include <string>
// #define _USE_MATH_DEFINES
#define LOG(x) cout<<x<<endl;
using namespace std;

class vec {
public:
	double x; double y; double z;

	double mod() {
		return sqrt(x * x + y * y + z * z);
	}
	void uguale(double val) {
		x = val; y = val; z = val;
	}
	void piu(double val) {
		x += val; y += val; z += val;
	}
	void per(double val) {
		x *= val; y *= val; z *= val;
	}
	void stampa() {
		cout << "x: " << x << "y: " << y << "z: " << z << "\n";
	}
};

#endif