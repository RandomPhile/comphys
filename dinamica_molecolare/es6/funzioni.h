#include "header.h"
struct vec {
    double x; double y; double z;

    double mod() {
        return sqrt(x * x + y * y + z * z);
    }
    void uguale(double val) {
        x = val; y = val; z = val;
    }
    void piu(double val) {
        x += val; y += val; z += val;
    }
    void per(double val) {
        x *= val; y *= val; z *= val;
    }
};
void vec_2D (vec *dr[N]) {
    for (int i = 0; i < N; ++i)
        dr[i] = new vec[N];
}
void stampa_stato(vec *r, vec *v, vec *a) {
    for (int i = 0; i < N; ++i) {
        cout << "i = " << i << endl;
        cout << "rx = " << r[i].x << endl;
        cout << "vx = " << v[i].x << endl;
        cout << "ax = " << a[i].x << "\n" << endl;

        cout << "ry = " << r[i].y << endl;
        cout << "vy = " << v[i].y << endl;
        cout << "ay = " << a[i].y << "\n" << endl;

        cout << "rz = " << r[i].z << endl;
        cout << "vz = " << v[i].z << endl;
        cout << "az = " << a[i].z << "\n" << endl;
        cout << "\n" << endl;
    }
}
void stampa_coord(vec *r) {
    for (int i = 0; i < N; ++i) {
        cout << "i = " << r[i].x << "\t" << r[i].y << "\t" << r[i].z << "\t" << endl;
    }
    LOG("\n");
}
void distr_gauss(vec *x, double sigma) {
    //Formule di Box-Muller
    for (int i = 0; i < N; ++i) {
        double x1, x2;
        x1 = rand() / ((double)RAND_MAX + 1.0);
        x2 = rand() / ((double)RAND_MAX + 1.0);

        x[i].x = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
        x[i].y = sigma * sqrt(-2 * log(1 - x2)) * sin(2 * M_PI * x1);

        x1 = rand() / ((double)RAND_MAX + 1.0);
        x2 = rand() / ((double)RAND_MAX + 1.0);

        x[i].z = sigma * sqrt(-2 * log(1 - x2)) * cos(2 * M_PI * x1);
    }
}
void v_cm_0(vec *v) {
    vec v_cm = {.x = 0, .y = 0, .z = 0};
    for (int i = 0; i < N; ++i) {
        v_cm.x += v[i].x;
        v_cm.y += v[i].y;
        v_cm.z += v[i].z;
    }
    v_cm.per(1 / N);
    for (int i = 0; i < N; ++i) {
        v[i].x -= v_cm.x;
        v[i].y -= v_cm.y;
        v[i].z -= v_cm.z;
    }
}
void plot_v_hist(vec *v) {
    //usare N grande
    double v_min = 1e100;
    double v_max = 0;
    double v_mod[N];
    for (int i = 0; i < N; ++i) {
        v_mod[i] = v[i].mod();
        if (v_mod[i] < v_min) {
            v_min = v_mod[i];
        } else if (v_mod[i] > v_max) {
            v_max = v_mod[i];
        }
    }
    int N_v = 1 + 3.322 * log(N); //Sturge's Rule
    double v_step = (v_max - v_min) / N_v;
    double bins;
    double f_v;
    for (int j = 0; j < N_v; ++j) {
        bins = v_min + v_step * j;
        f_v = 0;
        for (int i = 0; i < N; ++i) {
            if (v_mod[i] >= bins && v_mod[i] < (bins + v_step)) {
                f_v++;
            }
        }
        dati << bins << "\t" << f_v << endl;
    }
}