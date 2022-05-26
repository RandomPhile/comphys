#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include <armadillo>
#include <iomanip>
#include <algorithm>

#include "funz_gdr.h"
using namespace std;
using namespace arma;

int M=1;
int N=M*1000;

ofstream dati;

int main() {
	dati.open("dati.dat");

	double rho = 0.01;
	double L=cbrt(N/rho);
	cout<<L<<endl;
	mat r(N,3);

	int N_b=20;
	mat gdr(N_b, 2);//in una la gdr e nell'altra colonna la distanza dalla part centrale
	crea_reticolo(r, L);
	for (int i = 0; i < N; ++i){
		dati<<r(i,0)<<"\t"<<r(i,1)<<"\t"<<r(i,2)<<"\t"<<endl;
	}
	dati<<"\n\n";

	gdr_funz(r, L, rho, gdr, N_b);
	for (int i = 0; i < N_b; ++i){
		dati << gdr(i,1) << "\t" << gdr(i,0) << endl;
		//cout << gdr(i,1) << "\t" << gdr(i,0) << endl;
	}
	dati.close();
}

