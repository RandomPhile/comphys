#include "es11.h"

int main() {
    srand(1);//default seed = 1

    rowvec rho={0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
    rowvec passi_eq={400, 550, 1200, 1000, 2200, 2500, 2500, 2800, 3700, 3700, 3700, 3700, 3700, 3700, 3600, 2200};//tempi di equilibrazione con delta preso per avere circa 50%
    rowvec L(rho.size());
    for (int i = 0; i < rho.size(); ++i){
        L(i)=cbrt(N / rho(i));
    }

    ofstream gnuplot;
    gnuplot.open(gnuplot_path);
    gnuplot << caso_min << endl;
    gnuplot.close();

    // pressioni(rho, passi_eq, L);

    gdr_funz(L(caso_min), rho(caso_min));//da togliere il commento dopo aver già calcolato tutte le coordinate se caso_min!=-1

    return 0;
}

