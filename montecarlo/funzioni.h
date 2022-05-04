#ifndef funzioni_h
#define funzioni_h
#include "header.h"
#define _USE_MATH_DEFINES

double f(double x) {
    return exp(-x * x / 2);
}
double g(double x, double sigma, double mu) {
    return exp(-(x - mu) * (x - mu) / (2 * sigma * sigma));
}

void distr_unif(double *x, int len, double a, double b) {
    for (int i = 0; i < len; ++i) {
        x[i] = (b - a) * rand() / (RAND_MAX + 1.0) + a;
    }
}

void distr_gauss(double *x, int len, double sigma, double mu) {
    //Formule di Box-Muller
    double x1, x2;
    for (int i = 0; i < len; i += 2) {
        x1 = rand() / (RAND_MAX + 1.0);
        x2 = rand() / (RAND_MAX + 1.0);

        x[i]     = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1) + mu;
        x[i + 1] = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1) + mu;
    }
}

void distr_gauss_mezza(double *x, int len, double sigma, double mu) {
    //Formule di Box-Muller
    double x1, x2;
    for (int i = 0; i < len; i += 2) {
        x1 = rand() / (RAND_MAX + 1.0);
        x2 = rand() / (RAND_MAX + 1.0);

        x[i]     = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1) + mu;
        x[i + 1] = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1) + mu;
        if (x[i] < 3) {
            x[i] = x[i] + 2 * (3 - x[i]);
        }
        if (x[i+1] < 3) {
            x[i+1] = x[i+1] + 2 * (3 - x[i+1]);
        }
    }
    //se xi<3 "specchia" rispetto a 3
}
#endif