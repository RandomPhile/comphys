#include <iostream>
#include <cmath>
#include <array>
#include "funzioni.h"
using namespace std;

double max_of_function(double (*f)(double, double, int), int der, double a, double b, double param, int N_der);
double trapezi(double (*f)(double, double), double a, double b, int N, double param);
double simpson(double (*f)(double, double), double a, double b, int N, double param);

int main() {
    double A = 1, B = 1;
    double AB[] = {A, A, B, B};
    double alfa = 2, beta = 0.5;
    double param[] = {alfa, alfa, beta, beta};
    double I_t[4], I_s[4];
    int N_t[4], N_s[4];
    double f2_max[4], f4_max[4];


    double delta1=1e-4, delta2=1e-6;

    

    for (int i = 0; i < 4; ++i) {
        f2_max[i] = max_of_function(f_d[i], 2, 0, AB[i], param[i], 10);
        N_t[i] = ceil(sqrt(pow(AB[i], 3) * f2_max[i] / (12 * delta1)));

        f4_max[i] = max_of_function(f_d[i], 4, 0, AB[i], param[i], 10);
        N_s[i] = ceil(pow(pow(AB[i], 5) * f4_max[i] / (180 * delta2),0.25));
        printf("Max:%f,\t N_t:%d\n", f2_max[i],N_t[i]);
        printf("Max:%f,\t N_s:%d\n", f4_max[i],N_s[i]);
    }

    for (int i = 0; i < 4; ++i) {
        I_t[i] = trapezi(f[i], 0, AB[i], N_t[i], param[i]);
        I_s[i] = simpson(f[i], 0, AB[i], N_s[i], param[i]);
    }

    for (int i = 0; i < 2; ++i) {
        printf("I_t%d = %f\t\t", i + 1, I_t[i]);
        printf("I_s%d = %f\n", i + 1, I_s[i]);
    }
    printf("I_t%d = %f\t\t", 3, 1 / beta - I_t[2]);
    printf("I_s%d = %f\n", 3, 1 / beta - I_s[2]);
    printf("I_t%d = %f\t", 3, 2 / pow(beta, 3) - I_t[3]);
    printf("I_s%d = %f\n", 3, 2 / pow(beta, 3) - I_s[3]);

    return 0;
}

double max_of_function(double (*f)(double, double, int), int der, double a, double b, double param, int N_der) {
    float temp[N_der];
    double h_der;
    for (int n = 0; n < N_der; n++) {
        h_der = a + n * b / (1.0 * N_der);
        temp[n] = abs(f(h_der, param, der));
    }
    return *max_element(temp, temp + N_der);
}

double trapezi(double (*f)(double, double), double a, double b, int N, double param) {
    double h = (b - a) / N;

    double I;
    I = (f(a, param) + f(b, param)) / 2;
    for (int i = 1; i < N; i++) {
        I += f(a + i * h, param);
    }
    I *= h;
    return I;
}

double simpson(double (*f)(double, double), double a, double b, int N, double param) {
    double h = (b - a) / N;
    double I;
    I = (f(a, param) + f(b, param));
    for (int n = 1; n < N / 2; n++) {
        I += 2 * f(a + 2 * n * h, param) + 4 * f(a + (2 * n - 1) * h, param);
    }
    //tolgo l'elemento N/2 - 1 della prima sommatoria
    I = (I - 2 * f(a + 2 * (N / 2 - 1) * h, param)) * h / 3;
    return I;
}
