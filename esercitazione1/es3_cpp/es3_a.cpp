#include <iostream>
#include <cmath>

using namespace std;

void lista_n_float(float x_min, float x_max, int n, float *x) {
    float x_step = (x_max-x_min)*(1.f/n);
    for (int i = 0; i < n; ++i){x[i] = x_min + x_step*i;}
    return;
}
void lista_n_double(double x_min, double x_max, int n, double *x) {
    double x_step = (x_max-x_min)*(1.0/n);
    for (int i = 0; i < n; ++i){x[i] = x_min + x_step*i;}
    return;
}

/*
scelgo c t.c. il discriminante sia sempre positivo:
b^2 - 4ac >= 0 
c <= b^2/4a
*/

int main() {
    int n = 11;
    
    //float
    float x1[n], x2[n];
    float a = 1.f, b[n], b_min = 1, b_max = 11;
    float c = (b_min*b_min/(4.f*a))/1.1f;
    lista_n_float(b_min,b_max,n,b);
    
    //double
    double x1_[n], x2_[n];
    double a_ = 1.0, b_[n], b_min_ = 0.1, b_max_ = 10;
    double c_ = (b_min_*b_min_/(4.0*a_))*0.9;
    lista_n_double(b_min_,b_max_,n,b_);

    cout << "# Soluzione x1 (index 0)" << endl;

    for (int i = 0; i < n; ++i) {
        //float
        if (b[i]*b[i] - 4.f*a*c < 10e-8) {
            printf("Discriminante negativo = %f\n", b[i]*b[i] - 4.f*a*c);
        } else {
            x1[i] = (-b[i] + sqrt(b[i]*b[i] - 4.f*a*c))/(2.0*a);
            x2[i] = (-b[i] - sqrt(b[i]*b[i] - 4.f*a*c))/(2.0*a);
        }

        //double
        if (b_[i]*b_[i] - 4.f*a_*c_ < 10e-16) {
            printf("Discriminante negativo = %f\n", b_[i]*b_[i] - 4.f*a_*c_);
        } else {
            x1_[i] = (-b_[i] + sqrt(b_[i]*b_[i] - 4.0*a_*c_))/(2.0*a_);
            x2_[i] = (-b_[i] - sqrt(b_[i]*b_[i] - 4.0*a_*c_))/(2.0*a_);
        }
        cout << b[i] << "\t" << x1[i] << "\t" << x1_[i] << endl;
    }
    
    cout << "\n\n# Soluzione x2 (index 1)" << endl;
    for (int i = 0; i <  n; ++i) {
        cout << b[i] << "\t" << x2[i] << "\t" << x2_[i] << endl;
    }
    return 0;
}
