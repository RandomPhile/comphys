#ifndef es11_h
#define es11_h

#include "header.h"
#include "funzioni.h"

using namespace std;
using namespace arma;

extern int M;
extern int N;
extern int numero_accettati;
extern int numero_proposti;

extern ofstream dati;

double Calcola_Delta(double L, rowvec rho, int caso){
    double Delta;
    if(caso>=8){
        Delta=L/(60*rho(caso)*rho(caso));//scelgo un delta che mi dia circa 50% di accettazione
        cout<<"Delta scalato sulla scatola: "<<Delta/L<<endl;
        return Delta;
    }
    else if(caso<8 && caso>4){
        Delta=L/(50*rho(caso));//scelgo un delta che mi dia circa 50% di accettazione
        cout<<"Delta scalato sulla scatola: "<<Delta/L<<endl;
        return Delta;
    }
    else{
        Delta=L/(70*rho(caso));//scelgo un delta che mi dia circa 50% di accettazione
        if(caso!=-1){
            cout<<"Delta scalato sulla scatola: "<<Delta/L<<endl;
            return Delta;
        }
        else{
            return NAN;
        }
    }
}
void plot_pressioni(int caso_min) {
    if(caso_min==-1){
        string comando;
        comando = "gnuplot";
        comando += " plot2.plt";
        LOG(comando);
        system(comando.c_str());
    }
    else{
        string comando;
        comando = "gnuplot";
        comando += " plot.plt";
        LOG(comando);
        system(comando.c_str());
    }
}
double mod(rowvec &r){//calcolo il modulo della posizione relativa delle particelle i e j
    double mod= sqrt(pow1(r(0),2)+pow1(r(1),2)+pow1(r(2),2));
    return mod;
}
void crea_reticolo(mat &r, double L) {// passo la matrice per riferimento 
    int n = cbrt(N / M);
    double L_cella = L / cbrt(N / M);
    int cont = 0;//contatore particella
    mat b = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    if (M == 2) {b(1,2) = 0.5;}
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                rowvec R = {L_cella*(double)i, L_cella*(double)j, L_cella*(double)k};
                for (int l = 0; l < M; ++l) {
                    r.row(cont) = R + b.row(l) * L_cella;
                    cont++;
                }
            }
        }
    }
}
void gdr_plot(){
    //ora faccio il plot
    string comando;

    comando = "gnuplot";
    comando += " plot_gdr.plt";
    LOG(comando);
    system(comando.c_str());
}
void gdr_funz(mat &r, double L, double rho, int N_b){//penso funzionante
    cout<<"Ora calcolo la g(r), abbi ancora un po' di pazienza"<<endl;
    
    ofstream gdr_file;
    gdr_file.open("out/gdr_file.dat");
    
    mat gdr(N_b, 2);//in una la gdr e nell'altra colonna la distanza dalla part centrale        
    double delta_r=L/((double)N_b*2);
    rowvec dr(3);
    double dr_mod;

    for (int i = 0; i < N; ++i){//ciclo sulle particelle centrali
        for (int k = 0; k < N_b; ++k){//ciclo sui bin

            double R=delta_r*k+delta_r/2;//definisco il raggio medio del volumetto sferico
            double dV=4*M_PI*R*R*delta_r+M_PI/3*pow1(delta_r,3);//definisco il volumetto sferico
            double freq=0;//numero di particelle in un volumetto
    
            for (int j = 0; j < N; ++j){//ciclo sulle particelle non centrali
                
                for (int k1 = 0; k1 < 3; ++k1){
                    dr(k1) = r(i,k1) - r(j,k1);//trovo distanza relativa part i,j
                    dr(k1) -= L * rint(dr(k1)/L);//sposto in [-L/2,+L/2]
                }
                dr_mod=mod(dr);
                if(dr_mod<=R+delta_r/2 && dr_mod>R-delta_r/2){//se nel volumetto allora aumento la freq 
                    freq++;
                }
            }
            gdr(k,0)+=1/rho*(freq/dV);//def gdr
            gdr(k,1)=R;
        }
    }

    gdr.col(0)/=N;//medio sulle part
    for (int i = 0; i < N_b; ++i){
        gdr_file << (gdr(i,1)/L) * cbrt(N / M)<< "\t" << gdr(i,0) << endl;//in x c'è il raggio scalato sulla distanza iniziale dei primi vicini, la "cella"
    }
    gdr_plot();
}
void blocking_plot(){
    //ora faccio il plot
    string comando;

    comando = "gnuplot";
    comando += " plot_blocking.plt";
    LOG(comando);
    system(comando.c_str());
}
void blocking(int N_t){//crea e plotta il grafico del blocking
    ifstream dati;
    
    rowvec P(N_t);
    double variabile_inutile;

    ofstream blocking;
    blocking.open("blocking.dat");

    dati.open("dati.dat");
    
    for (int i = 0; i < N_t; ++i){
        dati >> variabile_inutile >> variabile_inutile >> variabile_inutile >> variabile_inutile >> P(i) >> variabile_inutile;//P(i) pressione istantanea
    }
    dati.close();

    for (int B = 1; B < N_t; B++){

        int N_B = floor(N_t / B);
        double P_mB=0;//pressione mediata sui blocchi
        double var_PB=0, P_media = 0;


        rowvec P_m(N_B, fill::zeros);//pressione media del blocco
        for (int i = 0; i < N_B; ++i){//ciclo sui blocchi per poi calcolare le medie
            for (int j = 0; j < B; ++j){
                P_m(i) += P(i * B + j) / B;//calcolo le medie sui blocchi
            }
            P_media += P(i) / N_t;
            P_mB += P_m(i) / (N_B);//calcolo la media complessiva come media sui blocchi
        }
        for (int i = 0; i < N_B; ++i){
            var_PB += (P_m(i) - P_mB) * (P_m(i) - P_mB) / (N_B);
        }
        // cout<<"Pressione media blocking "<<P_mB<<endl;
        blocking << B << "\t" << sqrt(var_PB / (N_B)) << endl;
    }

    blocking.close();
    blocking_plot();//faccio fare il plot
}

void blocking2(int N_t){//crea e plotta il grafico del blocking
    ifstream dati; ifstream risultati;
    
    rowvec P(N_t);
    double variabile_inutile;

    ofstream blocking;
    blocking.open("blocking.dat");

    dati.open("dati.dat");
    risultati.open("risultati.dat");
    
    for (int i = 0; i < N_t; ++i){
        dati >> variabile_inutile >> variabile_inutile >> variabile_inutile >> variabile_inutile >> P(i) >> variabile_inutile;//P(i) pressione istantanea
    }
    dati.close();
    risultati.close();


    for (int N_B = 0; N_B < N_t; ++N_B){
        int B = floor(N_t / N_B);
    }



    blocking.close();
    blocking_plot();//faccio fare il plot
}
void blocking(){
    for(int B = 1; B<t_max/10; B++){ // vario B, B particelle per blocco

    NB = floor(t_max/B); // NB è il numero di blocchi, floor approssima a intero
    // metto tutto a zero
    varpres=0;
    varenpot=0;
    pressionemedia=0;
    energiamedia=0;
    double mean[NB]; // vettore delle medie su ciascun blocco
    double meanen[NB]; // stessa cosa per pressione
    
    for (int jj = 1; jj <= NB; jj++){ // inizializzo a zero i vettori delle medie
        mean[jj]=0.;
        meanen[jj]=0.;
    }
    
    for (int jj = 0; jj < NB; jj++){ // ciclo sul numero di blocchi
        for (int kk=0; kk<B; kk++){ // ciclo sui B elementi del blocco
            mean[jj] += P[jj*B+kk]; // sommo tutte le pressioni istantanee del blocco jj esimo
            meanen[jj] += E[jj*B+kk]; // stessa cosa per le energie
        }
        mean[jj]=mean[jj]/B; // divido per B perché ho la media su B elementi, e questa è la media sul blocco jj esimo
        meanen[jj]=meanen[jj]/B; 
        pressionemedia += mean[jj]; // serve per la media totale
        energiamedia += meanen[jj];
    }
    pressionemedia = pressionemedia/NB; // questa è la media totale della pressione
    energiamedia=energiamedia/NB; //analogo per l'energia
    for (int i= 0; i < NB; i++){ // ciclo sui blocchi
        varpres += (mean[i]-pressionemedia)*(mean[i]-pressionemedia); 
        varenpot += (meanen[i]-energiamedia)*(meanen[i]-energiamedia);
    }
    varpres = varpres/NB; // è la varianza
    varenpot=varenpot/NB;
    
    // faccio l'errore finale 
    deltaP = sqrt(varpres/NB); // dalla formula
    printf( "%d\t %4.6e\n", B, deltaP); // per fare il grafico deltaP(B)
    deltaE=sqrt(varenpot/NB);
    //printf( "%d\t %4.6e\n", B, deltaE); // per fare il grafico deltaE(B)
}
}

#endif 
