/*
Nel file "integrali.h" abbiamo definito le 4 funzioni da integrare (più 2 test) con
le rispettive derivate, in modo che f1(x,a) = f1(x,a,0) ed f1''(x,a) = f1(x,a,2).
Nel file "funzioni.h" ci sono due funzioni per implementare l'integrazione numerica
col metodo dei trapezi e di Simpson, oltre ad una funzione per calcolare il massimo
di una funzione in un dato intervallo.
*/

#include <iostream>
#include <cmath>
#include "funzioni.h"
#include "integrali.h"
using namespace std;

int main() {
    double A = 1, B = 1, alpha = 2, beta = 0.5;

    double args[6][6] = {
        {0, A, alpha, 1e4, 0, 0},
        {0, A, alpha, 1e4, 0, 0},
        {0, B, beta , 1e4, 0, 0},
        {0, B, beta , 1e4, 0, 0},
        {-1,3, 0    , 1e4, 0, 0}, //funzione test retta
        {-2,2, 0    , 1e4, 0, 0}  //funzione test parabola
    };// a, b, param, N_fmax, N_t, N_s

    double delta1 = 1e-4, delta2 = 1e-6;
    double f_max, I_t, I_s;

    for (int i = 0; i < 6; ++i) {
        f_max = max_of_function(f[i], 2, args[i]); //massimo di f''
        args[i][4] = ceil(sqrt(pow(args[i][1] - args[i][0], 3.0) * f_max / (12.0 * delta1)));

        f_max = max_of_function(f[i], 4, args[i]); //massimo di f''''
        args[i][5] = ceil(pow(pow(args[i][1] - args[i][0], 5.0) * f_max / (180.0 * delta2), 1 / 4.0));

        switch (i) {
            //le funzioni test hanno derivata quarta nulla quindi N_s = 1
            //la retta ha anche derivata seconda nulla quindi N_t = 1
            case 4: args[i][4] = args[i][5] = 1; break;
            case 5: args[i][5] = 1; break;
        }

        I_t = trapezi(f[i], args[i]);//calcolo integrali
        I_s = simpson(f[i], args[i]);

        switch (i) {//utilizzo il risultato noto per l'integrale indefinito
        case 2:
            I_t = 1 / beta - I_t;
            I_s = 1 / beta - I_s;
            break;
        case 3:
            I_t = 2 / pow(beta, 3) - I_t;
            I_s = 2 / pow(beta, 3) - I_s;
            break;
        }

        //stampiamo i risultati
        printf("N_t%d = %3.0f\t", i + 1, args[i][4]);
        printf("N_s%d = %3.0f\t\t", i + 1, args[i][5]);
        printf("I_t%d = %9.6f\t", i + 1, I_t);
        printf("I_s%d = %9.6f\n", i + 1, I_s);
    }
    return 0;
}
