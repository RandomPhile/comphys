//
//  zeri.h
//  es2a
//
//  Created by Mattia Lupi on 16/03/22.
//

#ifndef zeri_h
#define zeri_h
extern std::ofstream dati;

double f_prime(double x);
double bisezione(double (*f)(double), double a, double b, double d){
    extern int contatore;
    double c;
    while (fabs(2 * (b - a) / (a + b)) > d) {
        c = (a + b) / 2;
        if (f(a)*f(c) < 0) {
            b = c;
        } else {
            a = c;
        }
        dati << contatore << '\t' << c << std::endl;
        contatore++;
    }
    return c;
}

double secante(double (*f)(double),double a, double b, double d) {
    double c;
    extern int contatore;
    while (fabs((b - a) / b) > d) {
        c = b - f(b) * (b - a) / (f(b) - f(a));
        a = b;
        b = c;
        dati << contatore << '\t' << c << std::endl;
        contatore++;
    }
    return c;
}

double NR(double (*f)(double), double a, int max) {
    extern int contatore;
    while (contatore <= max) {
        a = a - f(a)/f_prime(a);
        dati << contatore << '\t' << a << std::endl;
        contatore++;
    }
    return a;
}


#endif /* zeri_h */
