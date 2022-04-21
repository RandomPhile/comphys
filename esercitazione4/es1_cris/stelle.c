#include <stdio.h>
#include <math.h>
#include "eq_diff.h"

double f (double, double, double, double*); //identifico P con x, m con y e r con t
double g (double, double, double, double*);


int main () {

   //Inizializzo le condizioni iniziali e i vettori con i parametri
   double P_scale, R_scale, M_scale; //scale di dimensione
   P_scale = 150.1704; // MeV per fermi al cubo
   R_scale = 20.061455; //km
   M_scale = 13.6557594; //masse solari
   double P_dim, R_dim, M_dim; // P_dim = P_var * P_scale, ecc.
   
   double r = 1e-10; 
   double h = 1e-6;
   
   int N = 20; //numero di diverse pressioni iniziali
   double P_c[N][3]; //Ogni colonna è una stella, le righe sono le diverse pressioni iniziali
   double M[N][3];
   double R[N][3];
   double control[N][3]; //matrice di controllo dove calcolo M^(2 - gamma) * R^(3gamma - 4) e verifico che sia costante
   
   double fargs[] = {0.05, 5.0/3.0}; //i parametri di g e di f sono gli stessi, uso un'unica matrice. Ogni colonna è una stella
   
   //Queste sono condizioni che inizializzano in modo diverso P_c in base alla stella che stiamo studiando
    
      for (int j = 0; j < N; j++) { //NOTA: IN GENERALE, CONTROLLIAMO L'AFFIDABILITÀ DEI RISULTATI TRAMITE IL PARAMETRO CHE 
                                     //DOVREBBE ESSERE COSTANTE
         if ( j == 0 ) {
            P_c[j][0] = 0.00001;
            P_c[j][1] = 0.04;
            P_c[j][2] = 0.001;
         }
         else {
            P_c[j][0] =  2.514*P_c[j-1][0];
            P_c[j][1] =  1.484*P_c[j-1][1];
            P_c[j][2] =  1.543*P_c[j-1][2];
         }
      }
   
   /*printf("Stella 1: da %.5f a %.5f\n", P_c[0][0], P_c[N-1][0]);
   printf("Stella 2: da %.5f a %.5f\n", P_c[0][1], P_c[N-1][1]);
   printf("Stella 3: da %.5f a %.5f\n", P_c[0][2], P_c[N-1][2]);*/ //controllo gli estremi delle pressioni
  
   
   double P_var = P_c[0][0]; //la prima risoluzione per la prima stella la faccio "a mano" perché voglio i grafici P(r) e m(r)
   double m_var = 0;
   
   
   double *P, *m;
   P = &P_var;
   m = &m_var;
   
   FILE *st1;
   st1 = fopen("stella1.txt", "w");
   fprintf(st1, "%.5f\t%.5f\t%.5f\n", r*R_scale, P_var*P_scale, m_var*M_scale);
   
   
   
   while ( P_var > 1e-9 ) { 
   
      runge_kutta(r, h, P, m, f, fargs, g, fargs);
      r += h;    
      fprintf(st1, "%.5f\t%.5f\t%.5f\n", r*R_scale, P_var*P_scale, m_var*M_scale);
    }
    fclose(st1);
   
    M[0][0] = m_var; //massa e raggio finale della prima risoluzione della prima stella
    R[0][0] = r;
    control[0][0] = pow(M[0][0], 2 - fargs[1])*pow(R[0][0],3*fargs[1] - 4);
    
    //ora risolvo anche per gli altri valori di Pc, sempre per la stessa stella, e quindi stella colonna 0
    
    for (int j = 1; j < N ; j++ ) {
    
    r = 1e-10;
    P_var = P_c[j][0];
    m_var = 0;
    
       while (P_var > 1e-9) {
       
          runge_kutta(r, h, P, m, f, fargs, g, fargs);
          r += h;  
              
       }
    
    M[j][0] = m_var;
    R[j][0] = r;
    control[j][0] = pow(M[j][0], 2 - fargs[1])*pow(R[j][0],3*fargs[1] - 4);
      
    }
    
   //In questo modo ottengo i punti M(R) per la prima stella. Ora risolvo per le altre due
   
   for (int i = 1; i <= 2; i++) {
   
      if ( i == 1) {           //inizializzo gli argomenti in base a quale stella sto risolvendo, passaggio evitabile ma rende il codice più
                               //leggibile
         fargs[0] = 0.1;
         fargs[1] = 4.0/3.0;     
      }
      
      else {
         fargs[0] = 0.01;      
         fargs[1] = 2.54;
      }
      
      for (int j = 0; j < N; j++) {
      
         r = 1e-10;
         P_var = P_c[j][i];
         m_var = 0;
         
         while (P_var > 1e-9) {
         
            runge_kutta(r, h, P, m, f, fargs, g, fargs);
            r += h;
                    
         }
         
         M[j][i] = m_var;
         R[j][i] = r;
         control[j][i] = pow(M[j][i], 2 - fargs[1])*pow(R[j][i],3*fargs[1] - 4);      
      
      }   
   }
   
   FILE *mr;
   mr = fopen("massa_raggio.txt", "w");
   
   for (int i = 0; i <= 2; i++) { //ciclo sulle stelle
   
      for (int j = 0; j < N; j++) {
      
         if (j == (N - 1)) {
         
            fprintf(mr, "%.5f\t%.5f\t%.5f\n\n\n", R[j][i]*R_scale, M[j][i]*M_scale, control[j][i]); //separo i diversi blocchi
            
         }
         else {
         
            fprintf(mr, "%.5f\t%.5f\t%.5f\n", R[j][i]*R_scale, M[j][i]*M_scale, control[j][i]);
         
         }          
      }
   }   
      
      
      
   fclose(mr);
   
   return 0;
}








double f (double r, double P, double m, double *fargs ) {   //double fargs[] = {k, gam} 

   return - ( m / ( r * r ) ) * pow( P / ( fargs[0]*( fargs[1] - 1 ) ) , 1/fargs[1] );

}

double g (double r, double P, double m, double *gargs) {


   return r * r * pow( P / ( gargs[0] * ( gargs[1] - 1 ) )  , 1/gargs[1] );

}













