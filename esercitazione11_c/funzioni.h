#ifndef funzioni
#define funzioni

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define LOG(x) cout<<x<<endl;

void crea_reticolo(double *r, double L, int N, int M) {
	const int n = cbrt(N / M);
	const double L_cella = L / cbrt(N / M);

	int cont = 0;//contatore particella
	double b[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
	double R[3];

	if (M == 2) {b[1][2] = 0.5;}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			for (int k = 0; k < n; ++k) {
				R[0] = L_cella * (double)i;
				R[1] = L_cella * (double)j;
				R[2] = L_cella * (double)k;
				for (int l = 0; l < M; ++l) {
					r[cont + 0] = R[0] + b[l][0] * L_cella;
					r[cont + 1] = R[1] + b[l][1] * L_cella;
					r[cont + 2] = R[2] + b[l][2] * L_cella;
					cont += 3;
				}
			}
		}
	}
}

double V_LJ(double r, double L) { //potenziale
	double r6 = r * r * r * r * r * r;
	double L_mezzi = L / 2;
	double L_mezzi6 = L_mezzi * L_mezzi * L_mezzi * L_mezzi * L_mezzi * L_mezzi;
	return 4 * ((1 / (r6 * r6) - 1 / r6) - (1 / (L_mezzi6 * L_mezzi6) - 1 / L_mezzi6));
}
double dV_LJ(double r, double L) { //derivata potenziale
	double r6 = r * r * r * r * r * r;
	return -48 / (r6 * r6 * r) + 24 / (r6 * r);
}

#endif
