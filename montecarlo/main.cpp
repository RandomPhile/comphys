#include "header.h"

ofstream dati;

int main() {
    srand(1);
    int N = 5e5;
    double x[N], y[N];
    double a = 3, b = 1e5;
    double sigma = 0.5, mu = 3;
    distr_gauss_mezza(y, N, sigma, mu);

    distr_unif(x, N, a, b);

    double I, I2;
    for (int i = 0; i < N; ++i) {

        // LOG(x[i])
        I += f(x[i]);

        I2 += f(y[i]) / g(y[i], sigma, mu);
    }
    I *= (b - a) / N;

    double norm_g = (sqrt(2 * M_PI) * sigma)/2;
    LOG(norm_g)
    I2 *= norm_g / N;

    cout << "I_unif = " << I << endl;
    cout << "I_gaus = " << I2 << endl;
    cout << "Rif    = " << 0.003383692 <<endl;

}
