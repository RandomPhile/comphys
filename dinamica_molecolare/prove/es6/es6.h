#include "header.h"
void a_osc(vec *r, vec *a) {
    for (int i = 0; i < N; ++i) {//particella i=1..N
        a[i].x = -r[i].x;
        a[i].y = -r[i].y;
        a[i].z = -r[i].z;
    }
}
void vel_verlet(vec *r, vec *v, vec *a, double dt, double *K, double *V) {
    vec a_prev[N];
    for (int i = 0; i < N; ++i) {
        a_prev[i].x = a[i].x;
        a_prev[i].y = a[i].y;
        a_prev[i].z = a[i].z;

        r[i].x += v[i].x * dt + 0.5 * dt * dt * a[i].x;
        r[i].y += v[i].y * dt + 0.5 * dt * dt * a[i].y;
        r[i].z += v[i].z * dt + 0.5 * dt * dt * a[i].z;
    }
    a_osc(r, a);
    for (int i = 0; i < N; ++i) {
        v[i].x += 0.5 * dt * (a_prev[i].x + a[i].x);
        v[i].y += 0.5 * dt * (a_prev[i].y + a[i].y);
        v[i].z += 0.5 * dt * (a_prev[i].z + a[i].z);
    }
    
    *K = 0; *V = 0;
    double v_mod, r_mod;
    for (int i = 0; i < N; ++i) {
        v_mod = v[i].mod();
        r_mod = r[i].mod();

        *K += 0.5 * v_mod * v_mod;
        *V += 0.5 * r_mod * r_mod;
    }
}