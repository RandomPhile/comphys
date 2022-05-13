#include "funzioni.h"
#include <string>
/*** costanti globali ***/
const int M = 2; //1: CC, 2: BCC, 4: FCC¡
const int N = M * pow(4, 3);

const double T_req = 1.1;

int main() {
	/*** costanti simulazione ***/
	srand(2); //default seed = 1
	const double dt = 0.01; //passo temporale
	const double t1 = 10; //durata simulazione

	/*** costanti problema ***/
	coppia coppie[] = {//aggiornate con M=2,n=6,dt=0.01,t1=20 (12 minuti)
		{.rho = 0.01, .sigma = 1.04332},
		{.rho = 0.1, .sigma = 0.828379},
		{.rho = 0.2, .sigma = 0.670309},
		{.rho = 0.3, .sigma = 0.656069},
		{.rho = 0.4, .sigma = 0.762703},
		{.rho = 0.5, .sigma = 0.925801},
		{.rho = 0.6, .sigma = 1.12905},
		{.rho = 0.7, .sigma = 1.28666},
		{.rho = 0.8, .sigma = 1.41451},
		{.rho = 0.9, .sigma = 1.45276},
		{.rho = 1.0, .sigma = 1.44585},//
		{.rho = 1.1, .sigma = 1.39870},
		{.rho = 1.15, .sigma = 1.41261},
		{.rho = 1.2, .sigma = 1.45959}
	};

	string coord_path = "out/coordinate.xyz";
	string obs_path = "out/osservabili.dat";
	string p_path = "out/pressione.dat";
	string coord_g_path = "out/coordinate_gdr.xyz";
	string g_path = "out/gdr.dat";

	int caso = 5;//valore densità (e relativa sigma)
	
	calcolo_coordinate(coord_path, coppie[caso].rho, coppie[caso].sigma, dt, t1);
	calcolo_osservabili_da_file(coord_path, obs_path);

	int N_step    = 100;
	double pausa  = 0.01;
	// plot_coordinate(coord_path, N_step, pausa);
	plot_osservabili();

	// calcolo_pressioni(p_path, coppie, dt, t1, 0, 13);
	// plot_pressioni();

	int N_bins = 20;
	// calcolo_coordinate_per_gdr(coord_g_path, coppie[caso].rho, coppie[caso].sigma, dt, t1);
	// calcolo_gdr_da_file(coord_g_path, g_path, coppie[caso].rho, N_bins);
	// plot_gdr();

	return 0;
}