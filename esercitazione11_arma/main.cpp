#include "es11.h"

int main() {
    srand(1);//default seed = 1

    rowvec rho={0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
    rowvec passi_eq={400, 550, 1200, 1000, 2200, 2500, 2500, 2800, 3700, 3700, 3700, 3700, 3700, 3700, 3600, 2200};//tempi di equilibrazione con delta preso per avere circa 50%

    int caso_min = 0;//mettere -1 per avere P(rho)
    int caso_max;

    if (caso_min == -1) {
        caso_max = rho.size();
    } else {
        caso_max = caso_min + 1;
    }
    
    //###################################################
    mat r(N,3);
    cube dr(N,N,3);

    int reticolo   = log2(M);
    
    ofstream gnuplot;
    gnuplot.open(gnuplot_path);
    gnuplot << caso_min << endl;
    gnuplot.close();

    ofstream dati;
    dati.open(dati_path);

    ofstream coord;

    int start;
    if(caso_min==-1){
        start=caso_min+1;
    }
    else{
        start=caso_min;
    }
    for (int caso = start; caso < caso_max; ++caso) {

        N_t=passi_eq(caso)+1e4;//faccio fare 7000 passi dopo l'equilibrazione per avere dei risultati carini

        double L = cbrt(N / rho(caso));
        double Delta;
        double P = 0;

        Delta = calcola_delta(L, rho, caso);
        
        calcola_coord_oss(r, dr, L, passi_eq, caso, caso_min, Delta, &P);//da commentare una volta calcolate tutte le coord e osserv
        
        cout << "Densità = " << rho(caso) << endl;
        
        if (caso_min == -1) {
            dati << rho(caso) << "\t" << P << endl;
        }
        else{//calcolo della gdr
            // gdr_funz(coord_path, L, rho(caso));//da togliere il commento dopo aver già calcolato tutte le coordinate
        }
    }
    
    dati.close();

    return 0;
}

