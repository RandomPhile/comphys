#include "funzioni.h"

int main() {
	/*** costanti ***/
	const int M = 1; //1:CC, 2:BCC, 4:FCC
	const int N = M * pow(3, 3);
	const double T = 1.1;
	const int caso = 2; //tipo di sistema
	const int N_t  = 1e5; //t_eq + t_medie

	/*** variabili ***/
	srand(1); //default seed = 1
	struct {double rho, t_eq, delta;} sistema[] = {
		{0.1, 0, 1.91},
		{0.8, 0, 0.23},
		{1.2, 0, 0.11}
	}; double rho = sistema[caso].rho, t_eq = sistema[caso].t_eq, delta = sistema[caso].delta;
	if (t_eq >= N_t) {
		t_eq = N_t - 100;
	}
	double L = cbrt(N / rho);
	double *r = (double*)malloc(3 * N * sizeof(double)); //x1,y1,z1,x2,y2,z2

	/* variabili per MRT2 */
	int n; //molecola scelta in modo casuale ad ogni passo
	double r_n[3], dr[3], dr_n[3];
	double dr_mod, dr_n_mod;
	double A, dV;
	int numero_proposti = 0, numero_accettati = 0;

	double V, E, W, P, E_sum = 0, P_sum = 0;
	double var_E = 0, var_P = 0;
	const double K = 3 * N * T / 2;

	/* stampa informazioni sistema */
	fprintf(stderr, "M = %d\tN = %d\trho = %.2f\tL = %f\tL/2 = %f\n", M, N, rho, L, L / 2);

	/* variabili per blocking */
	int N_B; //numero blocchi
	FILE *dati_B = fopen("dati_blocking.dat", "w");

	for (int B = 1e1; B < 1e5; B *= 10) {
		N_B = N_t / B;
		fprintf(stderr, "\n%d blocchi da %d passi ciascuno\n", N_B, B);
		double *E_sum_B = (double*)malloc(N_B * sizeof(double));
		double *P_sum_B = (double*)malloc(N_B * sizeof(double));
		for (int b = 0; b < N_B; ++b) {
			E_sum_B[b] = 0; P_sum_B[b] = 0;
		}

		crea_reticolo(r, L, N, M);
		numero_proposti = 0, numero_accettati = 0;
		E_sum = 0, P_sum = 0; var_E = 0, var_P = 0;

		for (int t = 0; t < N_t; ++t) {

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

			E = K + V;
			P = rho * (1 + W / (3.0 * T));
			if (t >= t_eq) {
				E_sum += E;
				P_sum += P;
				var_E += E * E;
				var_P += P * P;
			}
			/* blocking */
			E_sum_B[t / B] += E;
			P_sum_B[t / B] += P;

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
		double P_avg = P_sum / (N_t - t_eq);
		double E_avg = E_sum / (N_t - t_eq);

		var_E = fabs((var_E / (N_t - t_eq)) - E_avg * E_avg) * (N_t - t_eq) / (N_t - t_eq - 1);
		var_P = fabs((var_P / (N_t - t_eq)) - P_avg * P_avg) * (N_t - t_eq) / (N_t - t_eq - 1);

		/* blocking */
		double var_E_B = 0, var_P_B = 0;
		for (int b = 0; b < N_B; ++b) {
			var_E_B += (E_sum_B[b] / B - E_avg) * (E_sum_B[b] / B - E_avg);
			var_P_B += (P_sum_B[b] / B - P_avg) * (P_sum_B[b] / B - P_avg);
		}
		var_E_B /= N_B;
		var_P_B /= N_B;

		fprintf(stderr, "P media finale = %f\t errore_P = %f\t", P_avg, sqrt(var_P / N_t));
		fprintf(stderr, "blocking: %f\t", sqrt(var_P_B / N_B));
		fprintf(stderr, "Tau_P = %f\n", B * var_P_B / var_P);

		fprintf(stderr, "E media finale = %f\t errore_E = %f\t", E_avg, sqrt(var_E / N_t));
		fprintf(stderr, "blocking: %f\t", sqrt(var_E_B / N_B));
		fprintf(stderr, "Tau_E = %f\n\n", B * var_E_B / var_E);

		fprintf(stderr, "Accettazione = %2.1f%%\tAccettati = %d\n", 100.0 * numero_accettati / (double)numero_proposti, numero_accettati);

		fprintf(dati_B, "%d\t%f\t%f\n", B, sqrt(var_P_B / N_B), sqrt(var_E_B / N_B));

		free(E_sum_B); free(P_sum_B);
	}
	fclose(dati_B); free(r);
	system("gnuplot plot_blocking.plt");
	return 0;
}