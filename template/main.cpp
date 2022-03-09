#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
ofstream dati;

int main() {
    dati.open("dati.dat");



    dati << "Prova dati" << endl;
    cout << "Prova cout" << endl;
    dati.close();
}
