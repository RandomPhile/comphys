#ifndef eq_diff_h
#define eq_diff_h
using namespace std;
extern int N; extern bool relativ;
extern ofstream dati, rect;
int eulero_exp(double t, double *x, double *y, double h,
               double (*f)(double, double, double, double*),
               double *fargs,
               double (*g)(double, double, double, double*),
               double *gargs) {
    double temp = *x;
    *x += h * f(t, *x, *y, fargs);
    *y += h * g(t, temp, *y, gargs);
    return 0;
}

int runge_kutta(double t, double *x, double *y, double h,
                double (*f)(double, double, double, double*),
                double *fargs,
                double (*g)(double, double, double, double*),
                double *gargs) {
    double k1 = h * f(t, *x, *y, fargs);
    double l1 = h * g(t, *x, *y, gargs);

    double k2 = h * f(t + h / 2, *x + k1 / 2, *y + l1 / 2, fargs);
    double l2 = h * g(t + h / 2, *x + k1 / 2, *y + l1 / 2, gargs);

    double k3 = h * f(t + h / 2, *x + k2 / 2, *y + l2 / 2, fargs);
    double l3 = h * g(t + h / 2, *x + k2 / 2, *y + l2 / 2, gargs);

    double k4 = h * f(t + h, *x + k3, *y + l3, fargs);
    double l4 = h * g(t + h, *x + k3, *y + l3, gargs);

    *x += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    *y += (l1 + 2 * l2 + 2 * l3 + l4) / 6;
    return 0;
}

double dP(double r, double P, double m, double *args) {
    double gam = args[0], K = args[1];
    double ro = pow(P / (K * (gam - 1)), 1.0 / gam);

    if (relativ) {
        double eps = ro + K * pow(ro, gam);
        return -(P + eps) * (m + r * r * r * P) / (r * r - 2 * m * r);
    } else {
        return -m * ro / (r * r);
    }
}

double dm(double r, double P, double m, double *args) {
    double gam = args[0], K = args[1];
    double ro = pow(P / (K * (gam - 1)), 1.0 / gam);

    if (relativ) {
        double eps = ro + K * pow(ro, gam);
        return r * r * eps;
    } else {
        return r * r * ro;
    }
}
void risolvi_stelle(double r0, double h, double *gam, double *K, double *Pc_max, double *Pc_min, double **M, double **Pc, double **R, bool funzione) {
    double m, P, r, cost;
    for (int i = 0; i < 3; ++i) {
        double args[] = {gam[i], K[i]};

        for (int j = 0; j < N; ++j) {
            Pc[i][j] = Pc_min[i] * pow( Pc_max[i] / Pc_min[i] , (double) j / (N - 1));

            P = Pc[i][j]; m = 0; r = r0;

            while (P > 0) {
                R[i][j] = r;
                M[i][j] = m;
                if (funzione == 0) {
                    runge_kutta(r, &P, &m, h, dP, args, dm, args);
                } else {
                    eulero_exp(r, &P, &m, h, dP, args, dm, args);
                }
                r += h;
            }
        }
    }
}
void stampa_valori(double **M, double **Pc, double **R, double R0, double M0, double *gam) {
    double cost;
    const char *colori[3] = {"red", "blue", "black"};
    const char *stelle[3] = {"{/Symbol G} = 5/3", "{/Symbol G} = 4/3", "{/Symbol G} = 2.54"};
    for (int i = 0; i < 3; ++i) {
        dati << "\n\n" << stelle[i] << endl;
        for (int j = 0; j < N; ++j) {
            cost = pow(M[i][j], 2 - gam[i]) * pow(R[i][j], 3 * gam[i] - 4);
            dati << Pc[i][j] << "\t" << R[i][j] << "\t" << M[i][j] << "\t" << cost << "\t" << R[i][j]*R0 << "\t" << M[i][j]*M0 << endl;
        }
    }
}

void stampa_valori_rel(double **M, double **Pc, double **R, double R0, double M0, double *gam) {
    rect.open("rect.dat");
    double deriv, cost;
    bool up = true;
    const char *colori[3] = {"red", "blue", "black"};
    const char *stelle[3] = {"{/Symbol G} = 5/3", "{/Symbol G} = 4/3", "{/Symbol G} = 2.54"};
    const char *comando = "\" fs transparent solid 0.15 border dashtype 2 lc rgb \"";
    for (int i = 0; i < 3; ++i) {
        // cout << "\n\nStella " << i + 1 << endl;
        dati << "\n\n" << stelle[i] << endl;
        for (int j = 0; j < N; ++j) {
            if (j < N - 1) {
                deriv = (M[i][j + 1] - M[i][j]) / (Pc[i][j + 1] - Pc[i][j]);
                deriv = deriv / abs(deriv);
                if (deriv > 0 and up == true) {
                    // cout << "inizio = " << Pc[i][j] << endl;
                    rect << "set object rect from " << Pc[i][j] << ",0 to ";
                    up = false;
                }
                if (deriv < 0 and up == false) {
                    // cout << "fine = " << Pc[i][j] << endl;
                    rect << Pc[i][j] << ",0.2 fc rgb \"" << colori[i] << comando << colori[i] << "\"" << endl;
                    up = true;
                }
            }

            cost = pow(M[i][j], 2 - gam[i]) * pow(R[i][j], 3 * gam[i] - 4);
            dati << Pc[i][j] << "\t" << R[i][j] << "\t" << M[i][j] << "\t" << cost << "\t" << R[i][j]*R0 << "\t" << M[i][j]*M0 << endl;
        }
        if (up == false) {
            // cout << "Fine = " << Pc[i][N - 1] << endl;
            rect << Pc[i][N - 1] << ",0.2 fc rgb \"" << colori[i] << comando << colori[i] << "\"" << endl;
            up = true;
        }
    }
    rect.close();
}

#endif /* eq_diff_h */