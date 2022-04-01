using namespace std;
extern double v;
extern int contatore;
extern ofstream dati;

double f(double x) {
    return 1 / tan(sqrt(v - x)) + sqrt(x / (v - x));
}

double f_prime(double x) {
    return 1 / (2 * pow(sin(sqrt(v - x)), 2) * sqrt(v - x)) + v / (2 * sqrt(x / (v - x)) * pow((x - v), 2));
}

double bisezione(double a, double b, double d) {
    double c;
    while (abs(2 * (b - a) / (a + b)) > d) {
        c = (a + b) / 2;
        if (f(a)*f(c) < 0) {
            b = c;
        } else {
            a = c;
        }
        dati << contatore << '\t' << c << endl;
        contatore++;
    }
    return c;
}

double secante(double a, double b, double d) {
    double c;
    while (abs((b - a) / b) > d) {
        c = b - f(b) * (b - a) / (f(b) - f(a));
        a = b;
        b = c;
        dati << contatore << '\t' << c << endl;
        contatore++;
    }
    return c;
}

double NR(double a, int max) {
    while (contatore <= max) {
        a = a - f(a) / f_prime(a);
        dati << contatore << '\t' << a << endl;
        contatore++;
    }
    return a;
}

double f1(double x, double a) {return pow(sin(a * x), 2);}
double f2(double x, double a) {return pow(x * sin(a * x), 2);}
double f3(double x, double b) {return exp(-b * x);}
double f4(double x, double b) {return pow(x, 2) * exp(-b * x);}

double simpson(double (*f)(double, double), double *fargs) {
    double a = fargs[0], b = fargs[1], param = fargs[2];
    int N = (int) fargs[3];

    double h = (b - a) / N;
    double I = 0;
    for (int n = 1; n <= (N / 2); n++) {
        I += f(a + 2 * (n-1) * h, param);
        I += 4.0 * f(a + (2 * n - 1) * h, param);
        I += f(a + 2 * n * h, param);
    }
    return I * (h / 3.0);
}