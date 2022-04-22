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
   
   double r = 1e-10; 
   double h = 1e-6;
   
   int N = 20; //numero di diverse pressioni iniziali
   
   double P_c[N][3]; //Ogni colonna è una stella, le righe sono le diverse pressioni iniziali
   double M[N][3];
   double R[N][3];
   double control[N][3]; //matrice di controllo dove calcolo M^(2 - gamma) * R^(3gamma - 4) e verifico che sia costante
   
   double P_var; //pressione e massa che uso come variabili
   double m_var;
   
   double *P, *m;
   P = &P_var;
   m = &m_var;
   
   double fargs[2]; //i parametri di g e di f sono gli stessi, uso un unico array.
   
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
            P_c[j][2] =  2*P_c[j-1][2]; //1.543
         }
      }
   
   /*printf("Stella 1: da %.5f a %.5f\n", P_c[0][0], P_c[N-1][0]);
   printf("Stella 2: da %.5f a %.5f\n", P_c[0][1], P_c[N-1][1]);
   printf("Stella 3: da %.5f a %.5f\n", P_c[0][2], P_c[N-1][2]);*/ //controllo gli estremi delle pressioni
   
   for (int i = 0; i <= 2; i++) { //con l'indice i ciclo sui diversi tipi di stelle, con j ciclo sulle condizioni iniziali
   
      if ( i == 0) { //risolvo per la prima stella; per la prima stella stampo anche i risultati in modo da eseguire un plot di
                     //P(r) e m(r)
                                                  
         fargs[0] = 0.05;
         fargs[1] = 5.0/3.0;
         
         char filename[15]; //stringa a cui sprintf assegnerà un nome diverso ad ogni ciclo
         FILE *st; //puntatore che uso per stampare i risultati su un file.txt
                  
         for (int j = 0; j < N ; j++ ) {
         
            r = 1e-10;
            P_var = P_c[j][0];
            m_var = 0;
            
            sprintf(filename, "stella%d.txt", j+1); //in modo da avere stella1.txt, ..., stellaN.txt
            st = fopen( filename, "w");
            fprintf(st, "%.5f\t%.5f\t%.5f\n", r*R_scale, P_var*P_scale, m_var*M_scale);
    
            while (P_var > 1e-9) {
       
               runge_kutta(r, h, P, m, f, fargs, g, fargs);
               r += h;  
               fprintf(st, "%.5f\t%.5f\t%.5f\n", r*R_scale, P_var*P_scale, m_var*M_scale);
              
            }
            fclose(st);
         
            M[j][i] = m_var;
            R[j][i] = r;
            control[j][i] = pow(M[j][i], 2 - fargs[1])*pow(R[j][i],3*fargs[1] - 4);
         
         }   
      
      }
   
      else if ( i == 1) {           //risolvo per la seconda stella
                                     
         fargs[0] = 0.1;
         fargs[1] = 4.0/3.0;    
         
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
      
      else  {         
                        //risolvo per la terza stella
         fargs[0] = 0.01;      
         fargs[1] = 2.54;
      
         
         for (int j = 0; j < N; j++) { 
      
      
            r = 1e-10;
            P_var = P_c[j][i];
            m_var = 0;
         
            while (P_var > 1e-7) {
         
               runge_kutta(r, h, P, m, f, fargs, g, fargs);
               r += h;
                    
            }
                 
            M[j][i] = m_var;
            R[j][i] = r;
            control[j][i] = pow(M[j][i], 2 - fargs[1])*pow(R[j][i],3*fargs[1] - 4);      
         }     
      
      }    
         
   }
   
   
   //Ho risolto tutto, ora stampo il file con massa totale in funzione di raggio totale, ossia gli elementi di M[j][i] e R[j][i]
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













