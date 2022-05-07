//
//  distribuzioni.h
//  
//
//  Created by Mattia Lupi on 07/05/22.
//

#ifndef distribuzioni_h
#define distribuzioni_h
void distr_unif(double *x, int L, double* args) {
    for (int i = 0; i < L; ++i) {
        x[i] = (args[1] - args[0]) * rand() / (RAND_MAX + 1.0) + args[0];
    }
}

void distr_gauss(double *x, int L, double *args) {
    //Formule di Box-Muller
    double x1, x2;
    for (int i = 0; i < L; i += 2) {
        x1 = rand() / (RAND_MAX + 1.0);
        x2 = rand() / (RAND_MAX + 1.0);

        x[i]     = args[2] * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1) + args[3];
        x[i + 1] = args[2] * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1) + args[3];
    }
}
void distr_gauss_mezza(double *x, int len, double *args) {
    //Formule di Box-Muller
    double x1, x2;
    for (int i = 0; i < len; i += 2) {
        x1 = rand() / (RAND_MAX + 1.0);
        x2 = rand() / (RAND_MAX + 1.0);

        x[i]     = args[2] * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1) + args[3];
        x[i + 1] = args[2] * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1) + args[3];
        if (x[i] < 3) {
            x[i] = x[i] + 2 * (3 - x[i]);
        }
        if (x[i+1] < 3) {
            x[i+1] = x[i+1] + 2 * (3 - x[i+1]);
        }
    }
    //se xi<3 "specchia" rispetto a 3
}


#endif /* distribuzioni_h */
