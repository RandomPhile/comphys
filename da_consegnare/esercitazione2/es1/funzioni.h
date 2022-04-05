using namespace std;
double max_of_function(double (*f)(double, double, int), int der, double *fargs) {
    int N = (int) fargs[3];
    float temp[N];
    double h;
    for (int i = 0; i < N; i++) {
        h = fargs[0] + i * fargs[1] / (1.0 * N);
        temp[i] = abs(f(h, fargs[2], der));
    }
    return *max_element(temp, temp + N);
}

double trapezi(double (*f)(double, double, int), double *fargs) {
    double a = fargs[0], b = fargs[1], param = fargs[2];
    int N = (int) fargs[4];

    double h = (b - a) / N;
    double I = (f(a, param, 0) + f(b, param, 0)) / 2.0;
    for (int n = 1; n < N; n++) {
        I += f(a + n * h, param, 0);
    }
    return I * h;
}

double simpson(double (*f)(double, double, int), double *fargs) {
    double a = fargs[0], b = fargs[1], param = fargs[2];
    if ((int)fargs[5] % 2 != 0) {fargs[5]++;}//se N è dispari aumento al prossimo pari
    int N = (int) fargs[5];

    double h = (b - a) / N;
    double I = 0;
    for (int n = 1; n <= (N / 2); n++) {
        I += f(a + 2 * (n-1) * h, param, 0);
        I += 4.0 * f(a + (2 * n - 1) * h, param, 0);
        I += f(a + 2 * n * h, param, 0);
    }
    return I * (h / 3.0);
}