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
int N=M*pow1(4,3);

ofstream dati;

void modifica_ret(mat &r,int righe, int colonne);

int main() {
	dati.open("dati.dat");

	double rho = 0.01;
	double L=cbrt(N/rho);
	cout<<L<<endl;
	mat r(N,3);

	int N_b=30;
	mat gdr(N_b, 2);//in una la gdr e nell'altra colonna la distanza dalla part centrale
	crea_reticolo(r, L);
	//modifica_ret(r,N,3);
	cout<<r<<endl;
	for (int i = 0; i < N; ++i){
		dati<<r(i,0)<<"\t"<<r(i,1)<<"\t"<<r(i,2)<<"\t"<<endl;
	}
	dati<<"\n\n";
	double media=0;
	gdr_funz(r, L, rho, gdr, N_b);
	for (int i = 0; i < N_b; ++i){
		dati << gdr(i,1) << "\t" << gdr(i,0) << endl;
		media+=gdr(i,0);
	}
	cout << media / N_b << endl;
	dati.close();
}

void modifica_ret(mat &r,int righe, int colonne){
	for (int i = 0; i < righe; ++i){
		for (int j = 0; j < colonne; ++j){
				r(i,j)=0;
			}	
	}
}