#include <iostream>
#include <cmath>
using namespace std;

typedef double (*lista_funzioni) (double x, double param);

double trapezi(double (*f)(double, double), double a, double b, int N, double param);
double simpson(double (*f)(double, double), double a, double b, int N, double param);

double f1(double x, double alfa) {
    return pow(sin(alfa * x), 2);
}
double f2(double x, double alfa) {
    return pow(x * sin(alfa * x), 2);
}
double f3(double x, double beta) {
    return exp(-beta * x);
}
double f4(double x, double beta) {
    return pow(x, 2) * exp(-beta * x);
}

int main() {
    lista_funzioni f[] = {f1, f2, f3, f4};
    double A = 1, B = 1;
    double AB[] = {A, A, B, B};
    double alfa = 2, beta = 0.5;
    double param[] = {alfa, alfa, beta, beta};
    double I_t[4], I_s[4];
    int N_t[4], N_s[4];

    for (int i = 0; i < 4; ++i) {
        N_t[i] = 500;
        N_s[i] = 4000000;
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
    for (int i = 1; i < N / 2; i++) {
        I = I + 2 * f(a + 2 * i * h, param) + 4 * f(a + (2 * i - 1) * h, param);
    }
    I = (I - 2 * f(b, param)) * h / 3;
    return I;
}
