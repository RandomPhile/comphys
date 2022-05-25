#include "costanti.h"
#include "funzioni.h"

int main() {
    double **a = new double*[N];
    for (int i = 0; i < N; ++i) {
        a[i] = new double[N];
    }

    double (*b)[N] = new double[N][N];

    double c[N][N];

    cout << typeid(a).name() << "\t\t" << typeid(a[0]).name() << "\t" << typeid(a[0][0]).name() << endl;
    cout << typeid(b).name() << "\t\t" << typeid(b[0]).name() << "\t" << typeid(b[0][0]).name() << endl;
    cout << typeid(c).name() << "\t\t" << typeid(c[0]).name() << "\t" << typeid(c[0][0]).name() << endl;

    return 0;
}
