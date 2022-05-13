#include "funzioni.h"

void crea_reticolo(vec *r, double L) {
	int n = cbrt(N / M);
	float a = L / cbrt(N / M);
	double b[][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
	if (M == 2) {b[1][2] = 0.5;}

	int cont = 0;//contatore particella
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			for (int k = 0; k < n; ++k) {
				vec R = {.x = i * a, .y = j * a, .z = k * a};
				for (int l = 0; l < M; ++l) {
					r[cont].x = R.x + b[l][0] * a;
					r[cont].y = R.y + b[l][1] * a;
					r[cont].z = R.z + b[l][2] * a;
					cont++;
				}
			}
		}
	}
}
void distr_gauss(vec *x, double sigma) {
	//Formule di Box-Muller
	for (int i = 0; i < N; ++i) {
		double x1, x2;
		x1 = rand() / ((double)RAND_MAX + 1.0);
		x2 = rand() / ((double)RAND_MAX + 1.0);

		x[i].x = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
		x[i].y = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1);

		x1 = rand() / ((double)RAND_MAX + 1.0);
		x2 = rand() / ((double)RAND_MAX + 1.0);

		x[i].z = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
	}
}
int input() {
	int input;
	cout << "\n"
	     "1. calcola coordinate per un valore di densità\n"
	     "2. plotta coordinate da file\n"
	     "3. calcola osservabili per un valore di densità\n"
	     "4. plotta osservabili da file\n";
	cout << "step: "; cin >> input;
	return input;
}
void vec_2D (vec *dr[], int dim) {
	for (int i = 0; i < dim; ++i)
		dr[i] = new vec[dim];
}
void aggiorna_a(vec *r, vec *a, double L) {
	vec dr;
	for (int i = 0; i < N; ++i) {
		a[i].uguale(0);

		for (int j = 0; j < N; ++j) {
			if (j != i) {
				dr.x = r[i].x - r[j].x;
				dr.y = r[i].y - r[j].y;
				dr.z = r[i].z - r[j].z;

				dr.x -= L * rint(dr.x / L);
				dr.y -= L * rint(dr.y / L);
				dr.z -= L * rint(dr.z / L);

				double dr_mod = dr.mod();

				if (dr_mod < L / 2) {
					double cost = -24 * (pow(1 / dr_mod, 8) - 2 * pow(1 / dr_mod, 14));

					a[i].x += dr.x * cost;
					a[i].y += dr.y * cost;
					a[i].z += dr.z * cost;
				}
			}
		}
	}
}
void stampa_coord(vec *r, vec *v, ofstream &file) {
	/* struttura file .xyz per VMD
	N
	nome molecola
	atomo1 x y z
	atomo2 x y z
	...
	atomoN x y z
	*/
	for (int i = 0; i < N; ++i) {
		file << "P" << i << "\t" << r[i].x << "\t" << r[i].y << "\t" << r[i].z << "\t" << v[i].x << "\t" << v[i].y << "\t" << v[i].z << "\n";
	}
}
void v_cm_0(vec *v) {
	vec v_cm = {.x = 0, .y = 0, .z = 0};
	for (int i = 0; i < N; ++i) {
		v_cm.x += v[i].x;
		v_cm.y += v[i].y;
		v_cm.z += v[i].z;
	}
	v_cm.per(1 / N);
	for (int i = 0; i < N; ++i) {
		v[i].x -= v_cm.x;
		v[i].y -= v_cm.y;
		v[i].z -= v_cm.z;
	}
}
void vel_verlet(vec *r, vec *v, vec *a, double dt, double L) {
	vec a_prev[N];
	for (int i = 0; i < N; ++i) {
		//POSIZIONI
		r[i].x += v[i].x * dt + 0.5 * dt * dt * a[i].x;
		r[i].y += v[i].y * dt + 0.5 * dt * dt * a[i].y;
		r[i].z += v[i].z * dt + 0.5 * dt * dt * a[i].z;

		//sposto dentro alla scatola
		// r[i].x -= L * floor(r[i].x / L);
		// r[i].y -= L * floor(r[i].y / L);
		// r[i].z -= L * floor(r[i].z / L);

		while (r[i].x > L) {r[i].x -= L;}
		while (r[i].y > L) {r[i].y -= L;}
		while (r[i].z > L) {r[i].z -= L;}

		while (r[i].x < 0) {r[i].x += L;}
		while (r[i].y < 0) {r[i].y += L;}
		while (r[i].z < 0) {r[i].z += L;}

		//copio le accelerazioni
		a_prev[i].x = a[i].x;
		a_prev[i].y = a[i].y;
		a_prev[i].z = a[i].z;
	}

	//ACCELERAZIONI
	aggiorna_a(r, a, L);

	for (int i = 0; i < N; ++i) {
		//VELOCITA
		v[i].x += 0.5 * dt * (a_prev[i].x + a[i].x);
		v[i].y += 0.5 * dt * (a_prev[i].y + a[i].y);
		v[i].z += 0.5 * dt * (a_prev[i].z + a[i].z);
	}
}
void calcolo_coordinate(string coord_path, double rho, double sigma, double dt, double t1) {
	double L = cbrt(N / rho);
	const int N_t = t1 / dt;
	vec r[N], v[N], a[N];//, *dr[N]; vec_2D(dr, N);
	ofstream coord;

	/*** inizializzo variabili ***/
	crea_reticolo(r, L);
	distr_gauss(v, sigma);
	aggiorna_a(r, a, L);

	coord.open(coord_path);
	coord << N << "\t" << N_t << "\t" << dt << "\t" << rho << "\t" << sigma << "\n";
	coord << "molecola\n";

	for (int i = 0; i < N_t; ++i) {
		stampa_coord(r, v, coord);
		vel_verlet(r, v, a, dt, L);
	}

	coord.close();
}
void plot_coordinate(string coord_path, int N_step, double pausa) {
	ifstream coord_in; string comando;
	int N_t; double rho, sigma, dt;

	coord_in.open(coord_path);
	coord_in >> N_t >> N_t >> dt >> rho >> sigma;
	coord_in.close();

	if (N_step > N_t) {N_step = N_t;}
	int skip = rint(N_t / N_step); if (skip == 0) {skip = 1;}
	double L = cbrt(N / rho);

	comando = "gnuplot -e N=";
	comando += to_string(N);
	comando += " -e L=";
	comando += to_string(L);
	comando += " -e N_step=";
	comando += to_string(N_step);
	comando += " -e dt=";
	comando += to_string(dt);
	comando += " -e skip=";
	comando += to_string(skip);
	comando += " -e pausa=";
	comando += to_string(pausa);
	comando += " animation.plt";
	LOG(comando);
	system(comando.c_str());
}
void calcolo_osservabili_da_file(string coord_path, string obs_path) {
	ifstream coord_in; ofstream obs;
	string line;
	int N_t; double rho, sigma, dt;

	coord_in.open(coord_path);
	coord_in >> N_t >> N_t >> dt >> rho >> sigma;
	coord_in >> line; //scarta la seconda riga

	double L = cbrt(N / rho);

	vec r[N], v, dr;
	double dr_mod, v_mod;

	obs.open(obs_path);
	double t = 0;
	double K_avg = 0, W_avg = 0, P, T;
	double K, V, E, W;
	for (int i_t = 0; i_t < N_t; ++i_t) {

		K = 0; V = 0; W = 0;

		for (int i = 0; i < N; ++i) {
			coord_in >> line; //scarta la prima colonna

			coord_in >> r[i].x >> r[i].y >> r[i].z;
			coord_in >> v.x >> v.y >> v.z;

			v_mod = v.mod();
			// obs << line <<endl;
			K += 0.5 * v_mod * v_mod;
		}

		for (int i = 0; i < N; ++i) {
			for (int j = i + 1; j < N; ++j) {
				dr.x = r[i].x - r[j].x;
				dr.y = r[i].y - r[j].y;
				dr.z = r[i].z - r[j].z;

				dr.x -= L * rint(dr.x / L);
				dr.y -= L * rint(dr.y / L);
				dr.z -= L * rint(dr.z / L);

				dr_mod = dr.mod();
				if (dr_mod < L / 2) {
					if (dr_mod < L / 2 && dr_mod != 0) {
						double VL_2 = 4 * (pow(2 / L, 12) - pow(2 / L, 6));
						//double VpL_2 = 24 * (pow(1 / dr_mod, 8) - 2 * pow(1 / dr_mod, 14));
						V += 4 * (pow(1 / dr_mod, 12) - pow(1 / dr_mod, 6)) - VL_2;
					}
					W += -24 * (pow(1 / dr_mod, 6) - 2 * pow(1 / dr_mod, 12));
				}
			}
		}
		W /= N;
		K_avg = K_avg + (K - K_avg) / (i_t + 1);//forse è +2
		W_avg = W_avg + (W - W_avg) / (i_t + 1);//forse è +2

		T = 2.0 * K_avg / (3.0 * N);
		P = (1 + W_avg / (3.0 * T_req));

		E = K + V;
		obs << t << "\t" << K << "\t" << V << "\t" << E << "\t" << T << "\t" << P << "\n";
		t += dt;
	}
	cout << "Rho = " << rho << "\t\tT = " << T << "\nSigma = " << sigma << "\t->\t" << sigma*sqrt(T_req / T) << "\n" << endl;
	coord_in.close();
	obs.close();
}
void plot_osservabili() {
	string comando;

	comando = "gnuplot";
	comando += " plot.plt";
	LOG(comando);
	system(comando.c_str());
}
double pressione(double rho, double sigma, double dt, double t1) {
	double L = cbrt(N / rho);
	const int N_t = t1 / dt;
	vec r[N], v[N], a[N], dr;//, *dr[N]; vec_2D(dr, N);
	double dr_mod;

	/*** inizializzo variabili ***/
	crea_reticolo(r, L);
	distr_gauss(v, sigma);
	aggiorna_a(r, a, L);

	double W_avg = 0, P;
	double W;

	double t = 0;
	for (int i_t = 0; i_t < N_t; ++i_t) {
		W = 0;
		vel_verlet(r, v, a, dt, L);

		for (int i = 0; i < N; ++i) {
			for (int j = i + 1; j < N; ++j) {
				dr.x = r[i].x - r[j].x;
				dr.y = r[i].y - r[j].y;
				dr.z = r[i].z - r[j].z;

				dr.x -= L * rint(dr.x / L);
				dr.y -= L * rint(dr.y / L);
				dr.z -= L * rint(dr.z / L);

				dr_mod = dr.mod();
				if (dr_mod < L / 2) {
					W += -24 * (pow(1 / dr_mod, 6) - 2 * pow(1 / dr_mod, 12));
				}
			}
		}
		W /= N;
		W_avg = W_avg + (W - W_avg) / (i_t + 1);//forse è +2
	}
	return (1 + W_avg / (3.0 * T_req));
}
void calcolo_pressioni(string p_path, coppia *coppie, double dt, double t1, int c_min, int c_max) {
	ofstream press;
	press.open(p_path);
	for (int caso = c_min; caso <= c_max; ++caso) {
		press << coppie[caso].rho << "\t" << pressione(coppie[caso].rho, coppie[caso].sigma, dt, t1) << "\n";
	}
	press.close();
}
void plot_pressioni() {
	string comando;

	comando = "gnuplot";
	comando += " plot_pressione.plt";
	LOG(comando);
	system(comando.c_str());
}
void calcolo_coordinate_per_gdr(string coord_g_path, double rho, double sigma, double dt, double t1) {
	ofstream coord_gdr;
	coord_gdr.open(coord_g_path);

	double L = cbrt(N / rho);
	const int N_t = t1 / dt;
	vec r[N], v[N], a[N], rp;
	int p = 0;//particella attorno cui calcolo g(r)

	/*** inizializzo variabili ***/
	crea_reticolo(r, L);
	distr_gauss(v, sigma);
	aggiorna_a(r, a, L);

	double t = 0;
	for (int i_t = 0; i_t < N_t; ++i_t) {
		vel_verlet(r, v, a, dt, L);
	}
	for (int i = 0; i < N; ++i) {
		rp.x = r[p].x - r[i].x;
		rp.y = r[p].y - r[i].y;
		rp.z = r[p].z - r[i].z;

		rp.x -= L * rint(rp.x / L);
		rp.y -= L * rint(rp.y / L);
		rp.z -= L * rint(rp.z / L);

		coord_gdr << rp.mod() << "\n";
	}
	coord_gdr.close();
}
void calcolo_gdr_da_file(string coord_g_path, string g_path, double rho, int N_bins) {
	ifstream coord_gdr;
	coord_gdr.open(coord_g_path);
	ofstream gdr;

	double L = cbrt(N / rho);

	double de_r = 0.5 * L / N_bins;
	cout << "de_r = " << de_r << endl;
	cout << "L/2 = " << L / 2 << endl;

	int freq[N_bins];
	for (int k = 0; k < N_bins; ++k) {
		freq[k] = 0;
	}
	double g, rp_mod, r_k;
	for (int i = 0; i < N; ++i) {
		coord_gdr >> rp_mod;
		for (int k = 0; k < N_bins; ++k) {
			r_k = (2 * k + 1) * de_r / 2;
			if (rp_mod > r_k - de_r / 2 && rp_mod <= r_k + de_r / 2) {
				freq[k]++;
				break;
			}
		}
	}
	coord_gdr.close();
	gdr.open(g_path);
	for (int k = 0; k < N_bins; ++k) {
		r_k = (2 * k + 1) * de_r / 2;
		g = freq[k] / ((pow(r_k + de_r / 2, 3) - pow(r_k - de_r / 2, 3)) * M_PI * 4.0 / 3.0);
		g /= rho;
		gdr << r_k << "\t" << freq[k] << "\t" << g << "\n";
	}
	gdr.close();
}
void plot_gdr() {
	string comando;

	comando = "gnuplot";
	comando += " plot_gdr.plt";
	LOG(comando);
	system(comando.c_str());
}