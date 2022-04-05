#include <stdio.h>
#include <math.h>

double f (double x); //è la funzione che integro fino a inf
double g (double x); //è la funzione che integro su un intervallo finito

double S_integra_f_inf(double a, double d); //integro con Simpson. a è l'estremo di partenza, d è la precisione
double T_integra_f_inf(double a, double d); //Stessa cosa, ma con trapezi
double S_integra_g (double a, double b, long int N); //integro g tra a e b con Simpson dividendo in N intervallini
double T_integra_g (double a, double b, long int N); //stessa cosa, con i trapezi

int S_counter = 0; //numero iterazioni Simpson
int T_counter = 0; //numero iterazioni trapezi
double alpha = 2.f;//alpha, beta, A, B, m ed n sono parametri che entrano nelle funzioni da integrare
double beta = 0.5;
double A = 1.f;
double B = 1.f;
int m = 0;
int n = 0;

int main () {

   double a_t, a_S, b_t, b_S, c_t, c_S, d_t, d_S;

   a_t = T_integra_g (0, A, 1e4);
   a_S = S_integra_g (0, A, 1e3);

   m = 2;

   b_t = T_integra_g (0, A, 1e4);
   b_S = S_integra_g (0, A, 1e3);

   c_t = T_integra_f_inf(B, 1e-7);
   c_S = S_integra_f_inf(B, 1e-6);

   n = 2;

   d_t = T_integra_f_inf(B, 1e-7);
   d_S = S_integra_f_inf(B, 1e-6);

   printf("   Trapezi\tSimpson\n");
   printf("a) %.5f\t%.5f\n", a_t, a_S);
   printf("b) %.5f\t%.5f\n", b_t, b_S);
   printf("c) %.5f\t%.5f\n", c_t, c_S);
   printf("d) %.5f\t%.5f\n", d_t, d_S);

   return 0;
}

double S_integra_f_inf(double a, double d) { //d è la precisione, a l'estremo iniziale

   double h = 0.1; //è il passo degli intervallini
   double b = 2 * h + a; //è l'estremo superiore dell'intervallino

   double I_0 = (h / 3) * (  f(a) + 4 * f((a + b) / 2)  + f(b)); //primo contributo a cui aggiungo tutti gli altri
   double I_n = I_0; //saranno i singoli contributi

   while ( fabs(I_n) > d ) {

      a = b;
      b = 2 * h + a;
      I_n = (h / 3) * ( f(a) + 4 * f((a + b) / 2)  + f(b));
      I_0 = I_0 + I_n;
      S_counter++;

   }

   return I_0;


}

double T_integra_f_inf(double a, double d) {

   double h = 0.01; //se integro con i trapezi, devo scegliere un passo più piccolo
   double b = a + h;

   double I_0 = (h / 2) * ( f(a) + f(b) ); //è quello in cui inglobo tutti i contributi
   double I_n = I_0; //saranno i singoli intervallini

   while ( fabs(I_n) > d ) {

      a = b;
      b = a + h;
      I_n = (h / 2) * ( f(a) + f(b) );
      I_0 = I_0 + I_n;
      T_counter++;

   }

   return I_0;
}

double S_integra_g (double a, double b, long int N) {

   double h = (b - a) / N;
   double I_0 = g(a) + g(b); //qui dentro ingolobo tutti i contributi

   for ( int i = 1 ; i <= (N / 2 - 1) ; i++ ) {

      I_0 = I_0 + 2 * g( a + 2 * i * h );

   }

   for ( int i = 1 ; i <= N / 2 ; i++) {

      I_0 = I_0 + 4 * g( a + (2 * i - 1) * h );

   }

   I_0 = (h / 3) * I_0;

   return I_0;

}

double T_integra_g (double a, double b, long int N) {

   double h = ( b - a ) / N;
   double I_0 = (h / 2) * ( g(a) + g(b) );

   for ( int k = 1 ; k <= ( N - 1) ; k++ ) {

      I_0 = I_0 + h * g( a + k * h );

   }

   return I_0;
}

double f (double x) {

   double z = pow(x, n) * exp(-beta * x);
   return z;
}

double g (double x) {

   double z = pow(x, m) * pow(sin(alpha * x), 2);
   return z;

}














