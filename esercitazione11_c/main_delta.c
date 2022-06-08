#include "funzioni.h"

int main() {
	/*** costanti ***/
	const int M = 4; //1:CC, 2:BCC, 4:FCC
	const int N = M * pow(5, 3);
	const double T = 1.1;
	const int N_t  = 1e5; //t_eq + t_medie

	/*** variabili ***/
	srand(1); //default seed = 1
	struct {double rho, t_eq, delta;} sistema[] = {//aggiornati con N=4*5^3, N_t=1e5
		{0.1, 0, 1.91},
		{0.8, 1e4, 0.23},
		{1.2, 2e4, 0.11}
	}; double rho, t_eq, delta;

	double L;
	double *r = (double*)malloc(3 * N * sizeof(double)); //x1,y1,z1,x2,y2,z2

	/* variabili per MRT2 */
	int n; //molecola scelta in modo casuale ad ogni passo
	double r_n[3], dr[3], dr_n[3];
	double dr_mod, dr_n_mod;
	double A, dV;
	int numero_proposti = 1, numero_accettati = 0; double rate = 1;

	int count;
	double incremento_delta = 0.01;
	for (int caso = 0; caso < 3; ++caso) {
		count = 0; rate = 1;
		rho = sistema[caso].rho, t_eq = sistema[caso].t_eq, delta = sistema[caso].delta;
		L = cbrt(N / rho);
		fprintf(stderr, "rho = %.2f\tdelta = %f\n", rho, delta);

		while (fabs(rate - 0.5) > 0.01 && count < 10) {
			numero_proposti = 0, numero_accettati = 0;
			crea_reticolo(r, L, N, M);
			for (int t = 1; t <= N_t; ++t) {
				/* metropolis */
				dV = 0;
				n = rint((N - 1) * (rand() / (RAND_MAX + 1.0))); //trovo la molecola che viene modificata da MTR2

				for (int c = 0; c < 3; ++c) { //x,y,z
					r_n[c] = r[3 * n + c] + delta * (rand() / (RAND_MAX + 1.0) - 0.5); //modifico le posizioni della particella n
				}

				for (int i = 0; i < 3 * N; i += 3) { //ogni altra particella (coppia contata una volta sola)
					if (i != 3 * n) {
						for (int c = 0; c < 3; ++c) { //x,y,z
							dr[c] = r[3 * n + c] - r[i + c];
							dr[c] -= L * rint(dr[c] / L);

							dr_n[c] = r_n[c] - r[i + c];
							dr_n[c] -= L * rint(dr_n[c] / L);
						}
						dr_mod   = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
						dr_n_mod = sqrt(dr_n[0] * dr_n[0] + dr_n[1] * dr_n[1] + dr_n[2] * dr_n[2]);

						dV += V_LJ(dr_n_mod, L) - V_LJ(dr_mod, L);
					}
				}

				A = exp(-dV / T);
				double num = rand() / (RAND_MAX + 1.0);
				numero_proposti++;
				if (A > num) {
					numero_accettati++;
					for (int c = 0; c < 3; ++c) { //x,y,z
						r[3 * n + c] = r_n[c];
					}
				}

			}
			rate = numero_accettati / (double)numero_proposti;
			// fprintf(stderr, "######nuovo delta = %f\tAccettati = %d (%2.1f%%)\n", delta, numero_accettati, 100.0 * rate);
			if (rate < 0.5) {
				delta -= incremento_delta;
			} else {
				delta += incremento_delta;
			}
			// rate = 0.5;//blocca il while
			count++;
		}
		rate = numero_accettati / (double)numero_proposti;
		if (0.011 < fabs(delta - sistema[caso].delta)) {
			fprintf(stderr, "######nuovo delta = %f\tAccettati = %d (%2.1f%%)\n", delta, numero_accettati, 100.0 * rate);
		}
	}
	free(r);
	return 0;
}