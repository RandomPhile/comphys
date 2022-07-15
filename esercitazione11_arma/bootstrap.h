void bootstrap(int N_t){//crea e plotta il grafico del blocking
    ifstream dati_blocking;
    
    rowvec P(N_t), E(N_t);

    ofstream bootstrap;
    bootstrap.open("out/bootstrap.dat");

    dati_blocking.open("out/dati_blocking.dat");
    
    for (int i = 0; i < N_t; ++i){
        dati_blocking >> P(i) >> E(i);//P(i) pressione istantanea, E(i) energia istantanea
    }
    dati_blocking.close();

    int N_B_prev=0;
    int N_boot = 200;

    for (int B = 300; B < N_t/3; B+=5){

        int N_B = floor(N_t / B);
        
        if(N_B!=N_B_prev){
            double P_media = 0, E_media = 0;

            cube PB(N_B, B, N_boot);
            cube EB(N_B, B, N_boot);

            mat PB_m(N_B, N_boot, fill::zeros);
            mat EB_m(N_B, N_boot, fill::zeros);

            rowvec PB_m_tot(N_B, fill::zeros);
            rowvec EB_m_tot(N_B, fill::zeros);

            for (int j = 0; j < N_B; ++j){//ciclo sui blocchi
                for (int i = 0; i < N_boot; ++i){//numero di bootstrap
                    for (int k = 0; k < B; ++k){//resampling dei punti per ogni blocco
                        int n = rint((rand() / (RAND_MAX + 1.)) * B);

                        PB(j, k, i) = P(n + j * B);//creo i blocchi resampled
                        EB(j, k, i) = E(n + j * B);

                        PB_m (j, i) +=  PB(j, k, i) / B;//calcolo le medie nei singoli blocchi 
                        EB_m (j, i) +=  EB(j, k, i) / B;
                    }

                    PB_m_tot(j) += PB_m (j, i) / N_boot;
                    EB_m_tot(j) += EB_m (j, i) / N_boot;
                }
                P_media += PB_m_tot(j) / N_B;
                E_media += EB_m_tot(j) / N_B;
            }

            double var_PB = 0;
            double var_EB = 0;

            for (int i = 0; i < N_B; ++i){
                var_PB += pow1(PB_m_tot(i) - P_media, 2) / N_B;
                var_EB += pow1(EB_m_tot(i) - E_media, 2) / N_B;
            }

            bootstrap << B << "\t" << sqrt(var_PB / N_B) << "\t" << sqrt(var_EB / N_B) << endl;

            N_B_prev=N_B;
        }
    }

    bootstrap.close();
}
