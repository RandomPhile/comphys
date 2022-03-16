//
//  funzioni.h
//  es2a
//
//  Created by Mattia Lupi on 16/03/22.
//

#ifndef funzioni_h
#define funzioni_h

extern double v;

typedef double (*lista_funzioni2) (double x, double param);
typedef double (*lista_funzioni3) (double x, double param, int der);

double f1(double x, double a) {return pow(sin(a * x), 2);}
double f2(double x, double a) {return pow(x * sin(a * x), 2);}
double f3(double x, double b) {return exp(-b * x);}
double f4(double x, double b) {return pow(x, 2) * exp(-b * x);}
lista_funzioni2 f[] = {f1, f2, f3, f4};

double f1_d(double x, double a, int der) {
    switch (der) {
    case 0:
        return f1(x, a);
    case 1:
        return a * sin(2 * a * x);
    case 2:
        return 2 * pow(a, 2) * cos(2 * a * x);
    case 3:
        return -4 * pow(a, 3) * sin(2 * a * x);
    case 4:
        return -8 * pow(a, 4) * cos(2 * a * x);
    }
    return 0;
}
double f2_d(double x, double a, int der) {
    switch (der) {
    case 0:
        return f2(x, a);
    case 1:
        return 2 * x * sin(a * x) * (sin(a * x) + a * x * cos(a * x));
    case 2:
        return 4 * a * x * sin(2 * a * x) + (2 * pow(a * x, 2) - 1) * cos(2 * a * x) + 1;
    case 3:
        return (6 * a - 4 * a * pow(a * x, 2) * sin(2 * a * x) + 12 * pow(a, 2) * x * cos(2 * a * x));
    case 4:
        return -8 * pow(a, 2) * (4 * a * x * sin(2 * a * x) + (pow(a * x, 2) - 3) * cos(2 * a * x));
    }
    return 0;
}
double f3_d(double x, double b, int der) {
    switch (der) {
    case 0:
        return f3(x, b);
    case 1:
        return -b * exp(-b * x);
    case 2:
        return pow(b, 2) * exp(-b * x);
    case 3:
        return -pow(b, 3) * exp(-b * x);
    case 4:
        return pow(b, 4) * exp(-b * x);
    }
    return 0;
}
double f4_d(double x, double b, int der) {
    switch (der) {
    case 0:
        return f4(x, b);
    case 1:
        return -x * (b * x - 2) * exp(-b * x);
    case 2:
        return (b * x * (b * x - 4) + 2) * exp(-b * x);
    case 3:
        return -b * (b * x * (b * x - 6) + 6) * exp(-b * x);
    case 4:
        return pow(b, 2) * (b * x - 6) * (b * x - 2) * exp(-b * x);
    }
    return 0;
}
lista_funzioni3 f_d[] = {f1_d, f2_d, f3_d, f4_d};

double max_of_function(double (*f)(double, double, int), int der, double a, double b, double param, int N_der) {
    float temp[N_der];
    double h_der;
    for (int n = 0; n < N_der; n++) {
        h_der = a + n * b / (1.0 * N_der);
        temp[n] = fabs(f(h_der, param, der));
    }
    return *std::max_element(temp, temp + N_der);
}
double f_prime(double x) {
    return 1 / (2 * pow(sin(sqrt(v - x)), 2) * sqrt(v - x)) + v / (2 * sqrt(x / (v - x)) * pow((x - v), 2));
}

#endif /* funzioni_h */
