#include <iostream>

/*
in pratica se scelgo n = 2m (un numero pari qualsiasi),
allora modificando un "altro_array", vado a modificare
il primo elemento dell' "array" originale.

okay, ho capito... il problema sono i <= invece di <
*/


int main() {
    int n = 2;
    float  array[n];
    double altro_array[n];

    for (int i = 0; i <= n; ++i) {
        array[i] = 0.2;
    }

    printf("%f\n\n", array[0]);

    for (int i = 0; i <= n; ++i) {
        altro_array[i] = 1;
        printf("%f\n", array[0]);
    }
    printf("\n%f\n", array[0]);
    
    return 0;
}
