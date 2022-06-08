#include "funzioni.h"

int main() {
	/*** costanti ***/
	const int M = 4; //1:CC, 2:BCC, 4:FCC
	const int N = M * pow(5, 3);
	const double T = 1.1;
	const int caso = 7; //tipo di sistema
	const int N_t  = 1e5; //t_eq + t_medie

	/*** variabili ***/
	srand(1); //default seed = 1
	struct {double rho, t_eq, delta;} sistema[] = {
		{0.1, 0, 1.91},
		{0.8, 1e4, 0.23},
		{1.2, 2e4, 0.11}
	}; double rho = sistema[caso].rho, t_eq = sistema[caso].t_eq, delta = sistema[caso].delta;
	if (t_eq >= N_t) {
		t_eq = N_t - 100;
	}
	double L = cbrt(N / rho);
	double *r = (double*)malloc(3 * N * sizeof(double)); //x1,y1,z1,x2,y2,z2
	crea_reticolo(r, L, N, M);

	/* variabili per MRT2 */
	int n; //molecola scelta in modo casuale ad ogni passo
	double r_n[3], dr[3], dr_n[3];
	double dr_mod, dr_n_mod;
	double A, dV;
	int numero_proposti = 0, numero_accettati = 0;

	/* gdr */
	FILE *dati_gdr = fopen("dati_gdr.dat", "w");
	const int N_bin = 30;
	const double de_r = 0.5 * 0.5 * L / N_bin;
	double r_k, dV_k;
	int freq[N_bin]; for (int i = 0; i < N_bin; ++i) {freq[i] = 0;}

	/* stampa informazioni sistema */
	fprintf(stderr, "M = %d\tN = %d\trho = %.2f\tL = %f\tL/2 = %f\tde_r = %f\n", M, N, rho, L, L / 2, de_r);

	for (int t = 1; t <= N_t; ++t) {

		/* metropolis */
		dV = 0;
		n = rint((N - 1) * (rand() / (RAND_MAX + 1.0))); //trovo la molecola che viene modificata da MTR2

		for (int c = 0; c < 3; ++c) { //x,y,z
			r_n[c] = r[3 * n + c] + delta * (rand() / (RAND_MAX + 1.0) - 0.5); //modifico le posizioni della particella n
		}

		for (int i = 0; i < 3 * N; i += 3) { //ogni altra particella
			if (i != 3 * n) {
				for (int c = 0; c < 3; ++c) { //x,y,z
					dr[c] = r[3 * n + c] - r[i + c];
					dr[c] -= L * rint(dr[c] / L);

					dr_n[c] = r_n[c] - r[i + c];
					dr_n[c] -= L * rint(dr_n[c] / L);
				}

				dr_mod   = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
				dr_n_mod = sqrt(dr_n[0] * dr_n[0] + dr_n[1] * dr_n[1] + dr_n[2] * dr_n[2]);

				//^^^ ho calcolato la distanza dr dalla particella n a tutte le altre
				for (int k = 0; k < N_bin; ++k) {
					r_k = (2 * k + 1) * de_r;
					if (dr_mod > r_k - de_r && dr_mod < r_k + de_r) {
						freq[k]++;
						break;
					}
				}

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

	for (int k = 0; k < N_bin; ++k) {
		freq[k] /= N_t;
		r_k = (2 * k + 1) * de_r;
		dV_k = (4.0 / 3) * M_PI * ((r_k + de_r) * (r_k + de_r) * (r_k + de_r) - (r_k - de_r) * (r_k - de_r) * (r_k - de_r));
		fprintf(dati_gdr, "%f\t%d\t%f\n", r_k, freq[k], freq[k] / (dV_k * rho));
	}
	fprintf(stderr, "Accettazione = %2.1f%%\tAccettati = %d\n", 100.0 * numero_accettati / (double)numero_proposti, numero_accettati);
	fclose(dati_gdr); free(r);
	system("gnuplot plot_gdr.plt");
	return 0;
}