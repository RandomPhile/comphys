#ifndef es11_h
#define es11_h

using namespace std;
using namespace arma;

extern int M;
extern int N;
extern int numero_accettati;
extern int numero_proposti;

extern ofstream dati;

#define LOG(x) cout<<x<<endl; 


void plot_reticolo() {
    string comando;
    comando = "gnuplot";
    comando += " plot_reticolo.plt";
    cout << comando << endl;
    system(comando.c_str());
    return;
}
void stampa_reticolo(mat &r){
    ofstream reticolo;
    reticolo.open("out/reticolo.dat");
    for (int i = 0; i < N; ++i){
        reticolo << r(i, 0) << "\t" << r(i,1) << "\t" << r(i,2) << endl;//0 posizione x, 1 posizione y, 2 posizione z
    }
    reticolo.close();
    // plot_reticolo();
    return;
}
double pow1(double base, int esp) {//funzione esponenziale creata per non usare pow
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}
double min(double a, double b){
    if(a<=b){
        return a;
    }
    else{
        return b;
    }
}
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
void plot_pressioni_energie(int caso_min) {
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
void crea_reticolo(mat &r, double L) {// passo la matrice per riferimento e crea il reticolo
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
void gdr_funz(double rho, int N_b){//penso funzionante
    mat r(N,3);

    ifstream reticolo;
    reticolo.open("out/reticolo.dat");
    for (int i = 0; i < N; ++i){//importa i dati dal reticolo
        reticolo >> r(i,0) >> r(i,1) >> r(i,2);
    }
    reticolo.close();

    double L = cbrt(N / rho);
    
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
    ifstream dati_blocking;
    
    rowvec P(N_t), E(N_t);

    ofstream blocking;
    blocking.open("out/blocking.dat");

    dati_blocking.open("out/dati_blocking.dat");
    
    for (int i = 0; i < N_t; ++i){
        dati_blocking >> P(i) >> E(i);//P(i) pressione istantanea, E(i) energia istantanea
    }
    dati_blocking.close();

    int N_B_prev=0;
    for (int B = 1; B < N_t/3; B+=10){

        int N_B = floor(N_t / B);//ceil dovuto all'uso di % N_t
        
        if(N_B!=N_B_prev){
            
            double P_mB=0;//pressione mediata sui blocchi
            double var_PB=0, P_media = 0;
            double E_mB=0;//energia mediata sui blocchi
            double var_EB=0, E_media = 0;
            
    
            rowvec P_m(N_B, fill::value(0));//pressione media del blocco
            rowvec E_m(N_B, fill::value(0));//energia media del blocco
            
            
            for (int i = 0; i < N_B; ++i){//ciclo sui blocchi per poi calcolare le medie
                for (int j = 0; j < B; ++j){
                    P_m(i) += P((i * B + j)) / B;//calcolo le medie sui blocchi 
                    E_m(i) += E((i * B + j)) / B;//calcolo le medie sui blocchi 
                }
                P_media += P(i) / N_t;
                P_mB  += P_m(i) / N_B;//calcolo la media complessiva come media sui blocchi
                E_media += E(i) / N_t;
                E_mB  += E_m(i) / N_B;//calcolo la media complessiva come media sui blocchi
            }
            for (int i = 0; i < N_B; ++i){
                var_PB += (P_m(i) - P_mB) * (P_m(i) - P_mB) / N_B;
                var_EB += (E_m(i) - E_mB) * (E_m(i) - E_mB) / N_B;
            }
            // cout<<"Pressione media blocking "<<P_mB<<endl;
            blocking << B << "\t" << sqrt(var_PB / N_B) << "\t" << sqrt(var_EB / N_B) << endl;

            N_B_prev=N_B;
        }
    }

    blocking.close();
}
void jackknife(int N_t){
    rowvec P(N_t), E(N_t);

    ifstream dati_blocking;
    ofstream jackknife;
    
    jackknife.open("out/jackknife.dat");
    dati_blocking.open("out/dati_blocking.dat");
    
    for (int i = 0; i < N_t; ++i){
        dati_blocking >> P(i) >> E(i);//P(i) pressione istantanea, E(i) energia istantanea
    }
    dati_blocking.close();

    int N_B_prev=1;//scarto il blocco singolo

    for (int B = (int)(N_t / 5e3) + 1; B < N_t/3; B+=5){//tolgo i blocchi troppo piccoli perche irrilevanti
        int N_B = floor(N_t / B);
        
        if(N_B!=N_B_prev){
            double P_media = 0, E_media = 0;

            rowvec P_m(N_B, fill::value(0));//pressione media del blocco
            rowvec E_m(N_B, fill::value(0));//energia media del blocco

            rowvec P_m_jack(N_B, fill::value(0));//pressione media del blocco
            rowvec E_m_jack(N_B, fill::value(0));//energia media del blocco

            double var_PB = 0;
            double var_EB = 0;

            //trovo le medie sui singoli blocchi per usarle dopo
            for (int i = 0; i < N_B; ++i){//ciclo sui blocchi per poi calcolare le medie
                for (int j = 0; j < B; ++j){//ciclo sui punti del blocco
                    P_m(i) += P((i * B + j)) / B;//calcolo le medie sui blocchi 
                    E_m(i) += E((i * B + j)) / B;//calcolo le medie sui blocchi 
                }
            }

            for (int i = 0; i < N_B; ++i){
                for (int j = 0; j < N_B; ++j){
                    if(j!=i){//scarto un blocco alla volta e faccio tutte le configurazioni
                        P_m_jack(i) += P_m(j) / (N_B - 1.0);
                        E_m_jack(i) += E_m(j) / (N_B - 1.0);
                    }
                }
                P_media += P_m_jack(i) / N_B;//trovo la media complessiva come media delle medie
                E_media += E_m_jack(i) / N_B;
            }

            for (int i = 0; i < N_B; ++i){
                var_PB += pow1(P_m_jack(i) - P_media, 2);
                var_EB += pow1(E_m_jack(i) - E_media, 2);
            }

            jackknife << B << "\t" << sqrt((N_B - 1) * var_PB / N_B) << "\t" << sqrt((N_B - 1) * var_EB / N_B) << endl;
            N_B_prev=N_B;
        }
    }

    jackknife.close();
}
#endif 
