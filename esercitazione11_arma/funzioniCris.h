#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.14159265

void crea_reticolo (double *x_var, double a, int n )  { //prende il vettore delle posizioni e il numero di particelle sugli spigoli

   int counter = 0; 
   
   for (int k = 0; k < n; k++) {
   
      for (int l = 0; l < n; l++) {
      
         for (int m = 0; m < n; m++) {
         
            x_var[3*counter] = a*k; //coordinata x della particella counter-esima
            x_var[3*counter + 1] = a*l; //stesso, y
            x_var[3*counter + 2] = a*m; //stesso, z
            
            counter++;
                    
         }     
       
      }
     
   }
     
}   


void distr_gauss(double *v_var, int N, double mu, double sigma) { //crea distr_gauss di media mu e dev. sigma

   double x_1, x_2;
   
   if ((3*N) % 2 == 0 ) { 
   
      for (int  j = 0; j < 3*N; j+=2){
      
         x_1 = rand() / (RAND_MAX + 1.0);
         x_2 = rand() / (RAND_MAX + 1.0);
         
         v_var[j] = sigma*sqrt( -2 * log( 1 - x_2 ) ) * cos( 2 * PI * x_1 ) + mu;
         v_var[j + 1] = sigma*sqrt( -2 * log( 1 - x_2 ) ) * sin( 2 * PI * x_1 ) + mu;
     
      }
   
   
   }
   
   else{
   
      for (int  j = 0; j < 3*N; j+=2){
      
         if (j == (3*N - 1)) {
         
            x_1 = rand() / (RAND_MAX + 1.0);
            v_var[j] = sigma*sqrt( -2 * log( 1 - x_2 ) ) * cos( 2 * PI * x_1 ) + mu;
            
         
         }
         
         else{
      
            x_1 = rand() / (RAND_MAX + 1.0);
            x_2 = rand() / (RAND_MAX + 1.0);
         
            v_var[j] = sigma*sqrt( -2 * log( 1 - x_2 ) ) * cos( 2 * PI * x_1 ) + mu;
            v_var[j + 1] = sigma*sqrt( -2 * log( 1 - x_2 ) ) * sin( 2 * PI * x_1 ) + mu;
         
         }
     
      }
      
   }

}


void setta_matr(double *vec, int dim, double val) {

   for (int j = 0; j < dim; j++ ) {
   
      vec[j] = val;
   
   }

}

void v_CM_0 (double *v_cm, double *v_var, int N) { //calcola la velocità del centro di massa e poi la sottrae a tutte quelle delle altre particelle


   for (int i = 0; i < N; i++) {
   
      v_cm[0] += v_var[3*i];
      v_cm[1] += v_var[3*i + 1];
      v_cm[2] += v_var[3*i + 2];
   
   }
   
   for (int j = 0; j < 3; j++) {
   
      v_cm[j] = v_cm[j] / N;
   }
   //Ora devo sottrarre a tutte le velocità questo valore
   
   for (int i = 0; i < N; i++) {
   
      v_var[3*i] = v_var[3*i] - v_cm[0];
      v_var[3*i + 1] = v_var[3*i + 1] - v_cm[1];
      v_var[3*i + 2] = v_var[3*i + 2] - v_cm[2];
   
   }



}

double V_LJ (double r, double L) {

   double r6 = r*r*r*r*r*r;
   double L_mezzi = L / 2.0;
   double L_mezzi_6 = L_mezzi * L_mezzi * L_mezzi * L_mezzi * L_mezzi * L_mezzi; 

   return 4 * (1.0 / (r6*r6) - 1.0 / r6 ) - 4 * (1.0 / (L_mezzi_6 * L_mezzi_6 ) - 1.0 / L_mezzi_6 );

}


double dV_LJ (double r) { //derivata potenziale

   double r6 = r*r*r*r*r*r;

   return -48.0 / (r6 * r6 * r) + 24.0 / (r6 * r);

}

void primi_vicini (double *x_var, double *dist, int N, double L, double r_c, int i) {


   for (int j = 0; j < N; j++) {
               
         dist[3*j] = ( x_var[3*i] - x_var[3*j] ) - L * rint( ( x_var[3*i] - x_var[3*j] ) / L );        
         dist[3*j + 1] = ( x_var[3*i + 1] - x_var[3*j + 1] ) - L * rint( ( x_var[3*i + 1] - x_var[3*j + 1] ) / L );        
         dist[3*j + 2] = ( x_var[3*i + 2] - x_var[3*j + 2] ) - L * rint( ( x_var[3*i + 2] - x_var[3*j + 2] ) / L );
                        
   } 

}


void calc_forze (double *x_var, double *dist, double *mod_dist, double *F, int N, double L, double r_c) {

   setta_matr(F, 3*N, 0);


  for (int i = 0; i < N; i++) {
           
      primi_vicini(x_var, dist, N, L, r_c, i); //fissata la particella i, calcolo le distanze con le altre con i primi vicini
      
         
      for (int j = 0; j < N; j++ ) {
      
         mod_dist[j] = sqrt( (dist[3*j] * dist[3*j]) + (dist[3*j + 1] * dist[3*j + 1]) + (dist[3*j + 2] * dist[3*j + 2]) );
      
      }
      //Adesso calcolo la forza che sente la particella i come somma delle forze lungo le direzioni
      
      
      
      for (int j = 0; j < N; j++) {      
   
         if (mod_dist[j] > 1e-4 && mod_dist[j] < r_c) { //evito una particella con se stessa
      
            F[3*i] -= dV_LJ (mod_dist[j]) * ( dist[3*j] / mod_dist[j] );
            F[3*i + 1] -= dV_LJ (mod_dist[j]) * ( dist[3*j + 1] / mod_dist[j] );
            F[3*i + 2] -= dV_LJ (mod_dist[j]) * ( dist[3*j + 2] / mod_dist[j] );
      
         }
   
   //Ho ottenuto le forze sulla particella i, ripeto per le altre   
      } 
    
   }
   
}


void blocking (double *P_vec, int N) { //N è il numero di punti, ossia di elementi in P_vec

   FILE *block;
   block = fopen("blocking.txt", "aw");
   
   int N_B_old = 0;
   

   
   for (int B = 1; B < N; B++) {
   
      if ((B % 10000) == 0) {
         printf("Blocking %d/%d\n", B, N);
      }
   
      int N_B = floor(N / (double) B ); //serve floor perché sennò rischio di leggere fuori da P_vec
      double P_med;
      double P_med_b[N_B];
      setta_matr(P_med_b, N_B, 0);
      double var_P = 0;
      double err_P;
      
      if (N_B_old != N_B) {
      
      for (int b = 0; b < N_B; b++){
         for (int i = 0; i < B; i++){
      
         P_med_b[b] += P_vec[b*B + i] / B; 
         
         }
         
      }   
      
      
      for (int j = 0; j < N; j++) {
      
         P_med += P_vec[j];
      
      }
      
      P_med /= N;
      
      for (int b = 0; b < N_B; b++) {
      
         var_P += (P_med_b[b] - P_med) * (P_med_b[b] - P_med) / N_B;
      
      }
      
      err_P = sqrt(var_P / N_B);
      
      fprintf(block, "%d\t%.5f\n", B, err_P );
      
      N_B_old = N_B;
      
      }
      }
   
    fclose(block);
   }
   
 































