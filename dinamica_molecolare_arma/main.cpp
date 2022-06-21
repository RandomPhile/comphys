#include "costanti.h"
#include "funzioni.h"

int main() {
	/*** costanti simulazione ***/
	srand(2); //default seed = 1
	const double dt = 0.01; //passo temporale
	const double t1 = 80; //durata simulazione

	/*** costanti problema ***/
	coppia coppie[] = {//eq di stato + gdr 
		{.rho = 0.01, .sigma = 1.04332, .t_eq = 15.2},//valore usato per relazione T,P,E
		{.rho = 0.1, .sigma = 0.828379, .t_eq = 10.1},
		{.rho = 0.2, .sigma = 0.670309, .t_eq = 8.3},
		{.rho = 0.3, .sigma = 0.656069, .t_eq = 4.6},
		{.rho = 0.4, .sigma = 0.762703, .t_eq = 2.7},
		{.rho = 0.5, .sigma = 0.925801, .t_eq = 2.2},
		{.rho = 0.6, .sigma = 1.12905, .t_eq = 2},
		{.rho = 0.7, .sigma = 1.28666, .t_eq = 1.6},
		{.rho = 0.8, .sigma = 1.41451, .t_eq = 1.2},//valore usato per relazione T,P,E
		{.rho = 0.9, .sigma = 1.45276, .t_eq = 0.94},
		{.rho = 1.0, .sigma = 1.44585, .t_eq = 0.55},
		{.rho = 1.1, .sigma = 1.39870, .t_eq = 0.51},
		{.rho = 1.15, .sigma = 1.41261, .t_eq = 0.51},
		{.rho = 1.2, .sigma = 1.45959, .t_eq = 0.62}//valore usato per relazione T,P,E
	};

	string coord_path = "out/coordinate.xyz";
	string obs_path = "out/osservabili.dat";
	string p_path = "out/pressione.dat";
	string coord_g_path = "out/coordinate_gdr.xyz";
	string g_path = "out/gdr.dat";

	int caso = 0;//valore densità (e relativa sigma) da utilizzare

	calcolo_coordinate(coord_path, coppie[caso].rho, coppie[caso].sigma, dt, t1);//calcola le coordinate, serve per calcolo osservabili e plot coordinate
	calcolo_osservabili_da_file(coord_path, obs_path, coppie[caso].t_eq);// calcola le osservabili, necessario per plot osservabili

	int N_step    = 100;
	double pausa  = 0.01;
	plot_coordinate(coord_path, N_step, pausa);//esegue l'animazione delle particelle secondo le coordinate calcolate
	plot_osservabili();//esegue i plot delle osservabili per il singolo caso di densità

	calcolo_pressioni(p_path, coppie, dt, t1, 0, 13, coppie[caso].t_eq);//esegue il calcolo delle pressioni per tutti i valori di densità dati
	plot_pressioni();//fa il plot delle pressioni in funzione della densità


	calcolo_coordinate_per_gdr(coord_g_path, coppie[caso].rho, coppie[caso].sigma, dt, t1);//calcola le coordinate per la g(r) per il caso di densità specificato 

	int N_bins = 35;
	calcolo_gdr_da_file(coord_g_path, g_path, coppie[caso].rho, N_bins);//calcola la g(r)
	plot_gdr();//esegue il plot della g(r)

	return 0;
}