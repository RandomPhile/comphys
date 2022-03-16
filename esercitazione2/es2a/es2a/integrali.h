#ifndef integrali_h
#define integrali_h
double simpson(double (*f)(double, double), double a, double b, int N, double param) {
    double h = (b - a) / N;
    //cout<<h<<endl;
    double I;
    I = (f(a, param) + f(b, param));
    for (int n = 1; n <= N / 2; n++) {
        I += 2 * f(a + 2 * n * h, param) + 4 * f(a + (2 * n - 1) * h, param);
    }
    //tolgo l'elemento N/2 della prima sommatoria
    I = (I - 2 * f(b, param)) * h / 3;
    return I;
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
#endif
