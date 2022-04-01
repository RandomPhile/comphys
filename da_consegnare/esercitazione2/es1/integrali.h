typedef double (*lista_funzioni) (double x, double param, int der);

double f1(double x, double a, int der) {
    switch (der) {
    case 0:
        return pow(sin(a * x), 2);
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
double f2(double x, double a, int der) {
    switch (der) {
    case 0:
        return pow(x * sin(a * x), 2);
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
double f3(double x, double b, int der) {
    switch (der) {
    case 0:
        return exp(-b * x);
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
double f4(double x, double b, int der) {
    switch (der) {
    case 0:
        return pow(x, 2) * exp(-b * x);
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

double f_test(double x, double b, int der) {
    switch (der) {
    case 0:
        return pow(x, 7) - pow(x, 6) + 2 * pow(x, 5) - 3 * pow(x, 4) + 5 * pow(x, 3) - 8 * pow(x, 2) + 13 * x - 21;
    case 1:
        return 7 * pow(x, 6) - 6 * pow(x, 5) + 10 * pow(x, 4) - 12 * pow(x, 3) + 15 * pow(x, 2) - 16 * x + 13;
    case 2:
        return 42 * pow(x, 5) - 30 * pow(x, 4) + 40 * pow(x, 3) - 36 * pow(x, 2) + 30 * x - 16;
    case 3:
        return 210 * pow(x, 4) - 120 * pow(x, 3) + 120 * pow(x, 2) - 72 * x + 30;
    case 4:
        return 840 * pow(x, 3) - 360 * pow(x, 2) + 240 * x - 72;
    }
    return 0;
}

lista_funzioni f[] = {f1, f2, f3, f4, f_test};