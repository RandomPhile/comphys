#include "costanti.h"
#include "funzioni.h"

int main() {
	/*** costanti simulazione ***/
	srand(2); //default seed = 1
	const double dt = 0.01; //passo temporale
	const double t1 = 40; //durata simulazione

	/*** costanti problema ***/
	coppia coppie[] = {//aggiornate con M=2,n=6,dt=0.01,t1=20 (12 minuti)
		{.rho = 0.01, .sigma = 1.04332, .t_eq = 0},//
		{.rho = 0.1, .sigma = 0.828379, .t_eq = 0},
		{.rho = 0.2, .sigma = 0.670309, .t_eq = 0},
		{.rho = 0.3, .sigma = 0.656069, .t_eq = 0},
		{.rho = 0.4, .sigma = 0.762703, .t_eq = 0},
		{.rho = 0.5, .sigma = 0.925801, .t_eq = 0},
		{.rho = 0.6, .sigma = 1.12905, .t_eq = 0},
		{.rho = 0.7, .sigma = 1.28666, .t_eq = 0},
		{.rho = 0.8, .sigma = 1.41451, .t_eq = 0},
		{.rho = 0.9, .sigma = 1.45276, .t_eq = 0},
		{.rho = 1.0, .sigma = 1.44585, .t_eq = 0},//
		{.rho = 1.1, .sigma = 1.39870, .t_eq = 0},
		{.rho = 1.15, .sigma = 1.41261, .t_eq = 0},
		{.rho = 1.2, .sigma = 1.45959, .t_eq = 0}
	};

	string coord_path = "out/coordinate.xyz";
	string obs_path = "out/osservabili.dat";
	string p_path = "out/pressione.dat";
	string coord_g_path = "out/coordinate_gdr.xyz";
	string g_path = "out/gdr.dat";

	int caso = 0;//valore densità (e relativa sigma)

	// calcolo_coordinate(coord_path, coppie[caso].rho, coppie[caso].sigma, dt, t1);
	// calcolo_osservabili_da_file(coord_path, obs_path);

	int N_step    = 100;
	double pausa  = 0.01;
	// plot_coordinate(coord_path, N_step, pausa);
	// plot_osservabili();

	calcolo_pressioni(p_path, coppie, dt, t1, 0, 13);
	plot_pressioni();


	// calcolo_coordinate_per_gdr(coord_g_path, coppie[caso].rho, coppie[caso].sigma, dt, t1);

	int N_bins = 30;
	// calcolo_gdr_da_file(coord_g_path, g_path, coppie[caso].rho, N_bins);
	// plot_gdr();

	return 0;
}