#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "eq_diff.h"

//per confrontare la convergenza dell'errore, partiamo da un passo grande, come 1e-2, 1e-3, e man mano divido per 2
double f (double, double, double, double*); //identifico P con x, m con y e r con t
double g (double, double, double, double*);

bool RK = 0; //mettere 1 se voglio la convergenza dell'errore con RK, viceversa mettere 0 se voglio la convergenza con Eulero

int main() {

   int N_E = 12; //numero di iterazioni, con RK bastano pochi punti (tre, quattro), con EULERO mettere almeno una decina
   double r, h_0, h; 
   int counter = 0; //contatore per le iterazioni
   int check, indice; //fattori che mi serve come condizione in un if
   double fargs[2]; //parametri di f e g
   double P_c = 47.367;
 
   
   r = 1e-10;
   h_0 = 1e-3; //step di partenza
   h = h_0;
   
   int dim0;                //so già quale sarà la dimensione dell'array di dati all'iterazione zero, chiaramente questo vale per h0 = 1e-3
   
   if (RK) {
      dim0 = 510;
   }
   else {
      dim0 = 505;
   }
   
   double vec_P_E[dim0][N_E + 1]; //matrice in cui salvo tutti i dati
     
   double P_var_E, m_var_E;
   double *P_E, *m_E;
   
   m_E = &m_var_E;
   P_E = &P_var_E;
   
  
   fargs[0] = 0.05; //inizializzo parametri e condizioni iniziali
   fargs[1] = 5.0/3.0;
                           
   FILE *st;
   
   char filename[20]; //stringhe a cui assegnerò nomi diversi ad ogni iterazione
 
   for (int n = 0; n <= N_E; n++) {
   
      P_var_E = P_c;
      m_var_E = 0;   
      r = 1e-10;
      
         
      sprintf(filename, "iterazione%d.txt", n);   
      st = fopen(filename, "w");
      
      counter = 0;
      
      
         if (n == 0) {
         
            while (P_var_E > 0) {
   
               fprintf(st, "%.5f\t%.5f\n", r, P_var_E);
               vec_P_E[counter][n] = P_var_E;
               counter++;  
               
               if (RK) {
                  runge_kutta(r, h, P_E, m_E, f, fargs, g, fargs);
               }
               else{     
                  eulero_exp(r, h, P_E, m_E, f, fargs, g, fargs);
               }
               r += h;
                                        
            }
               
         }
      
         else {  
         
            check = pow(2,n);
            
            while (P_var_E > 0) {
            
               if ( (counter % check ) == 0 && (counter/ (double) check ) < dim0 ) { //la prima volta non entra, poi dimezza l'indice
                                                                                     //per avere array della stessa dimensione nonstante dimezzo 
                                                                                     //il passo
      
                  fprintf(st, "%.5f\t%.5f\n", r, P_var_E);
                  indice = counter / (double) check;
                  vec_P_E[indice][n] = P_var_E;
     
               }
   
      
               counter++;       
               if (RK) {
                  runge_kutta(r, h, P_E, m_E, f, fargs, g, fargs);
               }
               else{     
                  eulero_exp(r, h, P_E, m_E, f, fargs, g, fargs);
               }
               r += h;
             
                  
            }
    
                   
         } 
         
      fclose(st);
         
      h /= 2;
   
   }
   //Quindi ora ho la matrice con dim0 righe e N_E + 1 colonne, e voglio confrontare tra le colonne, ossia le soluzioni via via più precise
     
   double err[dim0][N_E]; //facendo le differenze perdo una colonna; la colonna k di questa matrice è l'errore (relativo) della colonna k e k-1
                           //della matrice delle soluzioni
                           
   for (int n = 0; n < N_E; n++) {
   
      for (int j = 0; j < dim0; j++) {
      
      err[j][n] = fabs((vec_P_E[j][n + 1] - vec_P_E[j][n]) / vec_P_E[j][n]);
      
      }  
   }  
   
   //A questo punto calcolo un la correzione media rispetto alla soluzione precedente: devo fissare una colonna, sommare sulle righe e dividere per le righe
   
   double err_medi[N_E];
   double I; //variabile in cui accumulo le colonne
   
   for (int n = 0; n < N_E; n++) {
   
      I = 0;
   
      for(int j = 0; j < dim0; j++) {
      
         I += err[j][n];
      
      }
      
      err_medi[n] = I/dim0;
   
   }
   
   //A questo punto salvo gli errori medi in un file esterno per plottarli
   
   st = fopen("err.txt", "w");
   
   for(int n = 0; n < N_E; n++) {
   
   fprintf(st, "%.10f\t%.5f\n", h_0 / (pow(2,n+1)), err_medi[n]);
   
   }
   
   
   
     
   return 0;
}




double f (double r, double P, double m, double *fargs ) {   //double fargs[] = {k, gam} 
   
      return - ( m / ( r * r ) ) * pow( P / ( fargs[0]*( fargs[1] - 1 ) ) , 1/fargs[1] );
   
}

double g (double r, double P, double m, double *gargs) {


      return r * r * pow( P / ( gargs[0] * ( gargs[1] - 1 ) )  , 1/gargs[1] );
   

}









