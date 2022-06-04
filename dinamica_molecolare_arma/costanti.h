#ifndef COSTANTI_H
#define COSTANTI_H

#include <cmath>

#define ARMA_NO_DEBUG //non fa controlli e va piu veloce, da togliere solo se il codice funziona gia perfettamente


/*** costanti globali ***/
const int M = 1; //1: CC, 2: BCC, 4: FCC¡
const int N = M * pow(6, 3);

const double T_req = 1.1;

#endif
