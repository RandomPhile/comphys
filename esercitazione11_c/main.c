#include "funzioni.h"

int main() {
	/*** costanti ***/
	const int M = 1; //1:CC, 2:BCC, 4:FCC
	const int N = M * (const int)pow(3, 3);
	const double T = 1.1;
	const int caso = 7; //tipo di sistema
	const int N_t  = 1e5; //t_eq + t_medie

	/*** variabili ***/
	srand(1); //default seed = 1
	struct {double rho, t_eq, delta;} sistema[] = {
		{0.01, 0, 10},
		{0.1, 0, 1.91},
		{0.2, 1e4, 0.81},
		{0.3, 1e4, 0.55},
		{0.4, 1e4, 0.45},
		{0.5, 1e4, 0.38},
		{0.6, 1e4, 0.31},
		{0.7, 1e4, 0.27},
		{0.8, 1e4, 0.21},
		{0.9, 1e4, 0.18},
		{1.0, 1e4, 0.15},
		{1.1, 3e4, 0.12},
		{1.2, 2e4, 0.11}
	}; double rho = sistema[caso].rho, t_eq = sistema[caso].t_eq, delta = sistema[caso].delta;

	double L = cbrt(N / rho);
	double *r = (double*)malloc(3 * N * sizeof(double)); //x1,y1,z1,x2,y2,z2

	/* stampa informazioni sistema */
	fprintf(stderr, "M = %d\tN = %d\trho = %.2f\tL = %f\tL/2 = %f\n", M, N, rho, L, L / 2);

	crea_reticolo(r, L, N, M);

	// for (int i = 0; i < 3 * N; i += 3) fprintf(stderr, "%f\t%f\t%f\n", r[i + 0], r[i + 1], r[i + 2]);

	/* variabili per MRT2 */
	int n; //molecola scelta in modo casuale ad ogni passo
	double r_n[3], dr[3], dr_n[3];
	double dr_mod, dr_n_mod;
	double A, dV;
	int numero_proposti = 0, numero_accettati = 0;

	double V, E[N_t], W, P[N_t], E_sum = 0, P_sum = 0;
	const double K = 3 * N * T / 2;
	FILE *dati = fopen("dati.dat", "w");

	for (int t = 1; t <= N_t; ++t) {

		/* calcolo osservabili */
		V = 0; W = 0;
		for (int i = 0; i < 3 * N; i += 3) {
			for (int j = 0; j < i; j += 3) {
				for (int c = 0; c < 3; ++c) { //x,y,z
					dr[c] = r[i + c] - r[j + c];
					dr[c] -= L * rint(dr[c] / L);
				}
				dr_mod = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
				if (dr_mod <= L / 2) {
					V += V_LJ(dr_mod, L);
					W -= dr_mod * dV_LJ(dr_mod, L) / N;
				}
			}
		}

		E[t] = K + V;
		P[t] = rho * (1 + W / (3.0 * T));
		if (t >= t_eq) {
			E_sum += E[t];
			P_sum += P[t];
		}
		fprintf(dati, "%d\t%f\t%f\t%f\t%f\t%f\t\n", t, K, V, E[t], P[t], P_sum / (t - t_eq + 1));

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
		// fprintf(stderr, "A = %f\tnum = %f\n", A, num);
		if (A > num) {
			numero_accettati++;
			for (int c = 0; c < 3; ++c) { //x,y,z
				r[3 * n + c] = r_n[c];
			}
		}
	}

	// varianza campionaria
	double P_avg = P_sum / (N_t - t_eq + 1);
	double E_avg = E_sum / (N_t - t_eq + 1);
	double var_E = 0, var_P = 0;
	for (int t = t_eq; t < N_t; ++t) {
		var_E += (E[t] - E_avg) * (E[t] - E_avg);
		var_P += (P[t] - P_avg) * (P[t] - P_avg);
	}
	var_E /= N_t - t_eq - 1;
	var_P /= N_t - t_eq - 1;

	fprintf(stderr, "P media finale = %f\t sqrt(varianza) = %f\n", P_avg, sqrt(var_P));
	fprintf(stderr, "E media finale = %f\t sqrt(varianza) = %f\n", E_avg, sqrt(var_E));
	fprintf(stderr, "Accettazione = %2.1f%%\tAccettati = %d\n", 100.0 * numero_accettati / (double)numero_proposti, numero_accettati);
	fclose(dati); free(r);
	// system("gnuplot plot.plt");
	return 0;
}