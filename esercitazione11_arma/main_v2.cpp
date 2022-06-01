#include "costanti.h"
#include "funzioni.h"

int main() {
	/*** costanti simulazione ***/
	srand(2); //default seed = 1
	const double dt = 0.01; //passo temporale
	const double t1 = 40; //durata simulazione

	/*** costanti problema ***/
	rowvec rho={0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
    rowvec passi_eq={400, 550, 1200, 1000, 2200, 2500, 2500, 2800, 3700, 3700, 3700, 3700, 3700, 3700, 3600, 2200};//tempi di equilibrazione con delta preso per avere circa 50%

	string coord_path = "out/coordinate.xyz";
	string obs_path = "out/osservabili.dat";
	string p_path = "out/pressione.dat";
	string coord_g_path = "out/coordinate_gdr.xyz";
	string g_path = "out/gdr.dat";

	int caso = 0;//valore densità (e relativa sigma)

	double Delta;
    Delta=def_delta(rho,caso);//definisce il valore di delta in funzione del caso per avere accettazione circa 50%

	calcolo_coordinate(coord_path, rho(caso), passi_eq(caso), Delta);//sistemato per MRT2
	// calcolo_osservabili_da_file(coord_path, obs_path);//soon available for MRT2

	int N_step    = 100;
	double pausa  = 0.01;
	plot_coordinate(coord_path, N_step, pausa); //work in progress
	// plot_osservabili();//soon available for MRT2

	// calcolo_pressioni(p_path, coppie, dt, t1, 0, 13);//soon available for MRT2
	// plot_pressioni();//soon available for MRT2


	calcolo_coordinate_per_gdr(coord_g_path, rho(caso), coppie[caso].sigma, dt, t1);//soon available for MRT2

	int N_bins = 30;
	calcolo_gdr_da_file(coord_g_path, g_path, rho(caso), N_bins);//soon available for MRT2
	plot_gdr();//soon available for MRT2

	return 0;
}