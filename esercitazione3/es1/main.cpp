#include <iostream>
#include <fstream>
using namespace std;
ofstream dati;
//putatore NULL
double f(double t, double x, double y, double *fargs);
double g(double t, double x, double y, double *gargs);

int eulero_exp(double t, double *x, double *y, double h,
               double (*f)(double,double,double,double*),
               double *fargs,
               double (*g)(double,double,double,double*),
               double *gargs);

int main() {
    dati.open("dati.dat");
    double t ,t0 = 0, t1 = 25, h, x=1, y=0;
    int N = 100;
    h = (t1-t0)/( (double) N);
    printf("t\t\t\tx\t\t\ty\n");

    for (int n = 0; n < N; n++) {
        t = t0 + n*h;
        if (eulero_exp(t, &x, &y, h, f, NULL, g, NULL)){printf("ERRORE");}

        printf("%f\t%f\t%f\t\n", t,x,y);
        dati<<t<<"\t"<<x<<"\t"<<t<<endl;
    }
    dati.close();
    return 0;
}

double f(double t, double x, double y, double *fargs){
    return y;
};
double g(double t, double x, double y, double *gargs){
    return -x;
};

int eulero_exp(double t, double *x, double *y, double h,
               double (*f)(double,double,double,double*),
               double *fargs,
               double (*g)(double,double,double,double*),
               double *gargs){
    *x += h*f(t,*x,*y,NULL);
    *y += h*g(t,*x,*y,NULL);

    return 0;
}