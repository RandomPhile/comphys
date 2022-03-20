#include <stdio.h>
#include <math.h>

double bisezione (double a, double b, double d); //fa la bisezione tra a e b con precisione d
double secante (double a, double b, double d); //stesso discorso
double NR(double b, double d); //Newton-Raphson
double f(double x); //calcola la funzione f = cot(sqrt(v-e))+sqrt(e/(v-e)) nel punto in  e = x
double f_prime(double x); //calcola la derivata di f in x

double int_psi1 (double a, double b, long int N); //integra da 0 ad R con N intervallini
double int_psi2 (double a, double h, double d); //integra da a fino a +inf con intervallini grandi h

int counter_bisezione = 0; //conto quante iterazioni faccio per convergere alla soluzione
int counter_secante = 0;
int counter_NR = 0;

double v; //lo dichiaro qui perché deve essere visibile anche ad altre funzioni
double k;
double q;
int n = 2; //paramentro che entra nella funzione d'onda


int main() {

   double R = 1.93; //femtometri
   double V0 = 38.5; //MeV
   double Mpc2 = 938.272; //Mev
   double Mnc2 = 939.565; //Mev
   double muc2 = Mpc2*Mnc2/(Mpc2+Mnc2); //Mev
   double lam = pow(197.327,2)/(2*muc2*pow(1.93,2)); //Mev

   v = 2*V0*muc2*R*R/pow(197.327,2); //Parametro admensionale del problema che entra nella f
   printf("Il parametro v è pari a: %.5f\n", v); //Stampo v per controllo, ok circa 3.45 
   
   double e_bisezione = bisezione(0,3,1e-6); //zero di f con bisezione
   double E_bisezione = -lam *e_bisezione; //energia dello stato legato in MeV
   printf("Bisezione: e_star = %.5f, E = %.5f, %d iterazioni.\n", e_bisezione, E_bisezione, counter_bisezione); 
   
   double e_secante = secante(2.9,3,1e-6); //zero di f con secanti
   double E_secante = -lam*e_secante; //MeV
   printf("Secanti: e_star = %.5f, E = %.5f, %d iterazioni.\n", e_secante, E_secante, counter_secante);
   
   double e_NR = NR(3,1e-6);
   double E_NR = -lam*e_NR;
   printf("NR: e_star = %.5f, E = %.5f, %d iterazioni.\n", e_NR, E_NR, counter_NR);
   
   
   //Per calcolare la media di r^2 uso i risultati ottenuti con la bisezione, tanto sono uguali
   
   k = sqrt(2*muc2*(V0+E_bisezione)/pow(197.327,2)); //fermi alla meno uno
   q = sqrt(2*muc2*fabs(E_bisezione)/pow(197.327,2)); //fermi alla meno uno
   
   double N1 = int_psi1(0, R, 100);
   double N2 = exp(2*q*R)*pow(sin(k*R),2)*int_psi2(R, 0.1, 1e-5);
   
   n = 0;
   
   double D1 = int_psi1(0, R, 100);
   double D2 = exp(2*q*R)*pow(sin(k*R),2)*int_psi2(R, 0.1, 1e-5);
   
   double r_square = ( N1 + N2 ) / ( D1 + D2 );
   
   printf("Il valore d'aspettazione di r^2 sullo stato psi è: %.5f\n", r_square);
      
   return 0;
}



double bisezione (double a, double b, double d) { //applico bisezione e scrivo i risultati in un file txt

   double c = (a+b)/2;
   FILE *fbis;                                       //scrivo i passi e i risultati della bisezione in un file dati_bis.txt
   fbis = fopen("dati_bis.txt", "w");
   fprintf(fbis, "%d\t%.5f\n", counter_bisezione, c);
   fclose(fbis);

   while ( (2*fabs((b-a))/(b+a)) > d ) {
       
      if ( f(a)*f(c) < 0 ) {      
         b = c;
      }
      else {
         a = c;
      }
      
      c = (a+b)/2;
      counter_bisezione++;
      
      fbis = fopen("dati_bis.txt", "aw");
      fprintf(fbis, "%d\t%.5f\n", counter_bisezione, c);
      fclose(fbis);      
    }
    
    return c;
 }
 
 


 
double f(double x) { //sarebbe saggio mettere una condizione di verifica di positività delle radici ecc.
 
   double z = (1/tan(sqrt(v-x))) + sqrt(x/(v-x));
       
   return z;
}


double f_prime(double x) {

   double z = 1/(2*pow(sin(sqrt(v-x)),2)*sqrt(v-x)) + 0.5*sqrt((v-x)/x)*( v/( pow(v-x,2) ) );
   return z;

}


double psi1 (double r) {

   double x = pow(r,n)*pow(sin(k*r),2);
   return x;
}

double psi2 (double r) {

   double x = pow(r,n)*exp(-2*q*r);
   return x;

}
 
double secante (double a, double b, double d) {

   double c = b -  f(b) *  (b - a)/( f(b) - f(a) );
   
   FILE *fsec;
   fsec = fopen("dati_sec.txt","w");
   fprintf(fsec, "%d\t%.5f\n", counter_secante, c);
   fclose(fsec);
   
   while ( fabs((c-b)/c)  > d ) {
   
      a = b;
      b = c;
      
      c = b -  f(b) *  (b - a)/( f(b) - f(a) );
      counter_secante++;
      
      fsec = fopen("dati_sec.txt","aw");
      fprintf(fsec, "%d\t%.5f\n", counter_secante, c);
      fclose(fsec);
   
   }
   
   return c;


}

double NR (double b, double d) {

   double c = b - f(b)/f_prime(b);
   
   FILE *fNR;
   fNR = fopen("dati_NR.txt", "w");
   fprintf(fNR, "%d\t%.5f\n", counter_NR, c);
   fclose(fNR);
   
   
   while ( fabs(c - b)/c > d ) {
   
      b = c;
      
      c = b - f(b)/f_prime(b);
      counter_NR++;
      
      fNR = fopen("dati_NR.txt", "aw" );
      fprintf(fNR, "%d\t%.5f\n", counter_NR, c);
      fclose(fNR);
      
   }

   return c;

}

 
double int_psi1 (double a, double b, long int N) {

   double h = (b - a)/N;
   double I_0 = psi1(a) + psi1(b); //qui dentro ingolobo tutti i contributi
   
   for ( int i = 1 ; i <= (N/2 - 1) ; i++ ) {
   
      I_0 = I_0 + 2*psi1( a + 2*i*h );   
   
   }
   
   for ( int i = 1 ; i <= N/2 ; i++) {
   
      I_0 = I_0 + 4*psi1( a + (2*i - 1)*h );
   
   }
   
   I_0 = (h / 3) * I_0;
   
   return I_0;

}

double int_psi2(double a, double h , double d) { //d è la precisione, a l'estremo iniziale

   double b = 2*h + a;//è l'estremo superiore dell'intervallino
   
   double I_0 = (h / 3) * (  psi2(a) + 4*psi2((a+b)/2)  + psi2(b)); //primo contributo a cui aggiungo tutti gli altri
   double I_n = I_0; //saranno i singoli contributi
   
   while ( fabs(I_n) > d ) {
   
      a = b;
      b = 2*h + a;
      I_n = (h / 3) * ( psi2(a) + 4*psi2((a+b)/2)  + psi2(b));
      I_0 = I_0 + I_n;
    
   }
   
   return I_0;


}























     
      
      
      
      
      
      
      
      
