#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double V(double r, double *farg);

double ddr_V(double r, double *farg);

double gdr(double, double, double *, int , double *);

int main(void){
	
	// file
		FILE* fp_osservabili=fopen("osservabili.txt", "w+");
		FILE* fp_log=fopen("output.txt", "w+");
		//FILE* fp_medie=fopen("medie.txt", "w+");
		FILE *filep_gdr;
			char buffer_gdr[32];	
	
	// scale
		double mass=1.; //masse
		double eps=1.; //energie
		double sigma=1.; //lunghezze
		double T=sqrt(mass*sigma*sigma/eps); //tempo
	
	// parametri	
		double kbT0=1.1*eps; // fisso la T dell'NVT
		
		int M=4;
		int N=4*M*M*M; // qui fisso la N
		
		double dens=0.8 * sigma*sigma*sigma; //densità
		double L=pow(N/dens, 1./3); // qui fisso V (=L^3)
	
		double arg[5]={eps, sigma, L, dens, kbT0}; //vettore argomenti da passare alle funzioni
	
	// variabili
		double pos_arr[3*N];
		double dist_mtx[N*N], V_mtx[N*N], W_mtx[N*N];
		//double F_mtx[N*N*3]; //se volessi le forze
		double V_sistema=0., V_sistema_sum=0., V_sistema_mean=0.; //potenziale
		double W_sistema=0., W_sistema_sum=0., W_sistema_mean=0.; //contributo alla pressione
		double Press=0., Press_sum=0., Press_mean=0.; //pressione
	
	//cursori, contatori	
		int i, iii, j, jjj, k, l, lll, m, count; 
	
	// inizializzazione posizioni (FCC)
		double b[9]={0., 0.5, 0.5,
					 0.5, 0., 0.5,
					 0.5, 0.5, 0.};
		
		count=0;			 
		for(i=0;i<M;i++){
			for(j=0;j<M;j++){
				for(k=0;k<M;k++){
					
					// spigolo
					pos_arr[count]=i*L/M, count++;
					pos_arr[count]=j*L/M, count++;
					pos_arr[count]=k*L/M, count++;
					//printf("\n%d \t%lf \t%lf \t%lf", (count-3)/3, pos_arr[count-3], pos_arr[count-2], pos_arr[count-1]);
					
					// facce
					for(l=0;l<3;l++){
						lll=3*(l+1);
						for(m=0;m<3;m++){pos_arr[count]=pos_arr[count-lll]+b[3*l+m]*L/M, count++;}
						//printf("\n%d \t%lf \t%lf \t%lf", (count-3)/3, pos_arr[count-3], pos_arr[count-2], pos_arr[count-1]);	
					}
				}
			}
		}
		
	// calcolo distanze e potenziale
		double dist_comp[3]={0.,0.,0.};
		double ddrVij=0.;
		
		dist_mtx[0]=0., V_mtx[0]=0., W_mtx[0]=0.;
		for(i=1;i<N;i++){
		
			iii=3*i;
			for(j=0;j<i;j++){
				
				l=N*i+j, m=N*j+i; //puntano i coeff delle matrici
				dist_mtx[l]=0., dist_mtx[m]=0., V_mtx[l]=0., V_mtx[m]=0.; //inizializzo a 0
				
				jjj=3*j;
				for(k=0;k<3;k++){
					dist_comp[k] = pos_arr[iii+k] - pos_arr[jjj+k]; //la posizione di i nel sist ref di j
					dist_comp[k] = dist_comp[k] - L*rint(dist_comp[k]/L); //rint (primi vicini)
					//printf("\n pos_i=%lf \t pos_j=%lf \t dist_comp=%lf", pos_arr[iii+k], pos_arr[jjj+k], dist_comp[k]);
					
					dist_mtx[l] = dist_mtx[l] + dist_comp[k]*dist_comp[k];
				}
				dist_mtx[l]=sqrt(dist_mtx[l]);
				dist_mtx[m]=dist_mtx[l]; //carico la matrice
				//printf("\n dist_mtx[%d][%d]=%lf", i, j, dist_mtx[l]);
				
				if(dist_mtx[l]==0){
					printf("\nErrore! r[%d][%d]=0", i, j);
				}
				
				//calcolo il potenziale di interazione	
				V_mtx[l] = V(dist_mtx[l], arg);
				V_mtx[m] = V_mtx[l];
				//printf("\nV[%d][%d]=%lf", i, j, V_mtx[l]);
				V_sistema = V_sistema + V_mtx[l]; //sommo sulle coppie
				
				//calcolo le interazioni
				ddrVij = ddr_V(dist_mtx[l], arg);
				
				//e quindi il contributo alla pressione
				W_mtx[l] = W_mtx[l] - ddrVij * dist_mtx[l];
				W_mtx[m] = W_mtx[l];
				W_sistema = W_sistema + W_mtx[l]; //sommo sulle coppie
				
				/*
					//Se volessi le forze
					for(k=0;k<2;k++){
						F_mtx[l+k] = -ddrVij*dist_comp[k]/dist_mtx[l];
						F_mtx[m+k] = -F_mtx[l+k];
					}
				*/
			}
		}
		V_sistema=V_sistema/N; //divido per il numero di particelle
		W_sistema=W_sistema/N; //divido per il numero di particelle
		fprintf(fp_osservabili, "\n %d \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf", 0, V_sistema, V_sistema+1.5*kbT0, W_sistema, dens*W_sistema/3. + dens*kbT0);
		//printf("\n %d \t %10.8lf \t %10.8lf", 0, V_sistema, V_sistema+1.5*kbT0);
		
	// mando metropolis
		int part_num;
		int count_acc=0, count_den=0; // servono per conoscere il rate effettivo di accettazione
		
		double D=0.207; //pow(2, -10)
		double probs=0., randgenum=0.; // finestra
		
		double pos_curr[3]={0.,0.,0.}, pos_new[3]={0.,0.,0.};
		double dist_new[N], V_new[N], W_new[N];
		double Vpart_curr=0., Vpart_new=0.;
		double Wpart_curr=0., Wpart_new=0.;
		
		int bin_gdr=128;	//parametri per la g(r)
		double delta_r=L*0.5/bin_gdr;
		double r_k, g_k;
		
		srand(42);
		int passo, upbound=(7.6e+5), startmedie=1.4e+5, startpress=1.6e+5;
		for(passo=1;passo<upbound;passo++){
		
			// scelgo una particella
			part_num=passo%N;
			//printf("\nparticella %d", part_num);
			//printf("\n pos arr \t%lf \t%lf \t%lf", pos_arr[3*part_num], pos_arr[3*part_num+1], pos_arr[3*part_num+2]);
			
			// proposta
			for(k=0;k<3;k++){
				pos_curr[k]=pos_arr[3*part_num+k]; //salvo posizione corrente
				pos_new[k] = pos_curr[k] + D*( rand()/(RAND_MAX + 0.) - 0.5 ); //calcolo la mossa
			}
			//printf("\npos curr: \t%lf \t%lf \t%lf", pos_curr[0], pos_curr[1], pos_curr[2]);
			//printf("\npos new: \t%lf \t%lf \t%lf", pos_new[0], pos_new[1], pos_new[2]);
			
			// calcolo distanze e potenziali per proposta
			
			for(j=0;j<N;j++){
				
				jjj=3*j;
				dist_new[j]=0., V_new[j]=0.; //inizializzo a 0
				
				if(j!=part_num){
					
					for(k=0;k<3;k++){
						dist_comp[k] = pos_new[k] - pos_arr[jjj+k];
						dist_comp[k] = dist_comp[k] - L*rint(dist_comp[k]/L);
						
						dist_new[j] = dist_new[j] + dist_comp[k]*dist_comp[k];
					}
					dist_new[j]=sqrt(dist_new[j]);
					if(dist_new[j]==0){printf("\nErrore! dist_new[%d][%d]=0", part_num, j);}
					
					V_new[j]=V(dist_new[j], arg);
				}
			}
			
			Vpart_curr=0., Vpart_new=0.;
			for(j=0;j<N;j++){
				Vpart_curr = Vpart_curr + V_mtx[N*part_num+j];
				Vpart_new = Vpart_new + V_new[j];
			}
			fprintf(fp_log, "\n particella %d: \t Vc=%lf \t Vn=%lf \t Wc=%lf \t Wn=%lf", part_num, Vpart_curr, Vpart_new, Wpart_curr, Wpart_new);
			
			// calcolo probabilità di accettazione
			probs = exp( (Vpart_curr-Vpart_new)/kbT0 );
			//printf("\nprobs=%lf", probs);
			
			// lancio i dadi
			randgenum = rand()/(RAND_MAX + 1.);
			//printf("\nrand=%lf", randgenum);
			
			if(randgenum<probs){ //se probs è maggiore di 1 lo sarà in automatico di randgenum
				
				for(i=0;i<3;i++){pos_arr[3*part_num+i]=pos_new[i];} //aggiorno posizione
				//printf("\npos = %lf \t %lf \t %lf", pos_arr[part_num+0], pos_arr[part_num+1], pos_arr[part_num+2]);
				
				Wpart_curr=0., Wpart_new=0.;
				for(j=0;j<N;j++){
					jjj=3*j;
					l=N*part_num+j; m=N*j+part_num;
					W_new[j]=0.; //inizializzo a 0
					
					if(j!=part_num){
						ddrVij = ddr_V(dist_new[j], arg);
						W_new[j] = W_new[j] - ddrVij * dist_new[j];
					}
					Wpart_curr = Wpart_curr + W_mtx[l];
					Wpart_new = Wpart_new + W_new[j];
					
					dist_mtx[l]=dist_new[j], dist_mtx[m]=dist_mtx[l]; //aggiorno distanze
					
					V_mtx[l]=V_new[j], V_mtx[m]=V_mtx[l]; //aggiorno potenziale (mtx)
					
					W_mtx[l]=W_new[j], W_mtx[m]=W_mtx[l]; //aggiorno il contributo alla pressione
					//printf("\n dist[%d]=%lf \t V[%d]=%lf", j, dist_mtx[l], j, V_mtx[l]);
				}
				
				V_sistema = N*V_sistema - Vpart_curr + Vpart_new; //aggiorno potenziale
				V_sistema = V_sistema/N;
				W_sistema = N*W_sistema - Wpart_curr + Wpart_new; //aggiorno il contributo alla pressione
				W_sistema = W_sistema/N;
				fprintf(fp_log,"\nVsist=%10.8lf \t Wsist=%10.8lf\n", V_sistema, W_sistema);
				count_acc++;	//incremento il counter delle accettazioni
			}
			else{ 
					fprintf(fp_log,"\n mossa rifiutata: passo %d", passo);
					fprintf(fp_log,"\nVsist=%10.8lf \t Wsist=%10.8lf\n", V_sistema, W_sistema);
					count_den++;
				} //incremento il counter dei rifiuti
			fprintf(fp_osservabili, "\n %d \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf", passo, V_sistema, V_sistema+1.5*kbT0, W_sistema, dens*W_sistema/3. + dens*kbT0);
			//printf( "\n %d \t %10.8lf \t %10.8lf", passo, V_sistema, V_sistema+kbT0);
			
		/*------------------CALCOLO DELLE MEDIE----------------------------------------------------------------------------------------*/	
			
			if(passo>startmedie){ //per il calcolo delle medie
				V_sistema_sum = V_sistema_sum + V_sistema;
				W_sistema_sum = W_sistema_sum + W_sistema;
			}
			
		/*------------------CALCOLO DELLA G(R)----------------------------------------------------------------------------------------*/
			/*
			if(passo>startmedie && passo%N==3){
				sprintf(buffer_gdr, "gdr_%d.txt", passo);
				filep_gdr=fopen(buffer_gdr, "w+");
				
				for(k=0;k<bin_gdr;k++){
					r_k=(2*k+1.)*delta_r*0.5;
					g_k=gdr(r_k, delta_r, dist_mtx, N, arg);
					fprintf(filep_gdr, "\n %10.8lf \t %10.8lf", r_k, g_k);
				}
				fclose(filep_gdr);
			}
			*/	
		}
	
	/*---------------------------STAMPA DELLE MEDIE------------------------------------------------------------------------------------*/
		
		V_sistema_mean=V_sistema_sum/(passo-startmedie);
		Press_mean=dens*W_sistema_sum/(3.*(passo-startmedie)) + kbT0*dens;
		
		fprintf(fp_log,"\ncount_acc=%d \tcount_den=%d", count_acc, count_den);
		printf("\ncount_acc=%d \tcount_den=%d", count_acc, count_den);
		fprintf(fp_log,"\nV_mean=%10.8lf \tE_mean=%10.8lf \tW_mean=%10.8lf \tPress_mean=%10.8lf", 
					V_sistema_mean, V_sistema_mean+1.5*kbT0, W_sistema_mean, Press_mean);

		fclose(fp_osservabili);
		
	/*------------------CALCOLO DELLA VARIANZA------------------------------------------------------------------------------------*/
	
	//devo calcolare in primis i coefficienti di correlazione
		
		
		return 0;
}

double V(double r, double *farg){
	
	double eps=farg[0];
	double sigma=farg[1];
	double Lm=farg[2]*0.5;

	//printf("\n \n Chiamata a U con valore r=%lf", r);
	double dummy1 = (sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r);
	//printf("\n \t sig/r^6=%lf", dummy1);
	//double dummy2 = dummy1*dummy1;
	//printf("\n \t sig/r^12=%lf", dummy2);
	//printf("\n \t ritorna valore %lf \n", 4*eps * ( dummy2-dummy1 ));
	return 4.*eps * dummy1*( dummy1 - 1. );
	
}

double ddr_V(double r, double *farg){
	
	// memo: Harg[4]={eps, sigma, L, dens};
	double eps=farg[0];
	double sigma=farg[1];
	double Lm=farg[2]*0.5;
	
	if(r>Lm){return 0;}
	
	else{
		//printf("\n \n Chiamata a U con valore r=%lf", r);
	double dummy1 = (sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r);
	//printf("\n \t sig/r^6=%lf", dummy1);
	double dummy2 = dummy1*dummy1;
	//printf("\n \t sig/r^12=%lf", dummy2);
	//printf("\n \t ritorna valore %lf \n", 4*eps * ( dummy2-dummy1 ));
	return -4.*eps * ( 12.*dummy2 - 6.*dummy1 )/r;
	}
}

double gdr(double r, double delta_r, double *rij_gdr, int N, double *arg){
	
	//printf("\n r=%lf \t delta_r=%lf", r, delta_r);
	double rp, rm, Vol=0, Num=0;
	rp=r+delta_r*0.5, rm=r-delta_r*0.5;
	int i=0, j=0, l=0, count=0;
	for(i=1;i<N;i++){
		for(j=0;j<i;j++){
			
			l=N*i+j;
			if( rij_gdr[l]>rm && rij_gdr[l]<rp ){count++;}
		}	
	}
	
	/*
		questi conteggi tengono in considerazione tutte le N*(N-1)/2 interazioni (rij=rji in modulo)
		questo richiede un fattore di normalizzazione N/2 nel calcolo di g(r)
		-- caso limite per g(r0)=1 mi aspetterei di contare N/2 distanze con valore ~r0 --
		-- (se conto che j dista r0 da i poi non conto che i dista r0 da j) --
	*/
	Num=count*2./N; 
	//printf("\n \t %d", Num);
	
	Vol=4*M_PI*( rp*rp*rp - rm*rm*rm )/3;
	//printf("\n \t Vol=%lf", Vol);
	
	//printf("\n \t gdr=%lf", Num/(Vol*arg[3]) );
	return Num/(Vol*arg[3]); //arg[3] è dens
}
