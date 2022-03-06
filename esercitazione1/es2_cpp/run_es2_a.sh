#!/bin/bash
g++ es2_a.cpp -o es2_a.out
./es2_a.out > dati_a.dat
pkill -x gnuplot, gnuplot_qt
#!/usr/bin/gnuplot -persist
gnuplot -persist <<-EOFMarker

set termoption enhanced
set title "Esercizio 2a"
set xlabel "h"
set ylabel "\{/Symbol D}_h"
set grid
set logscale x
plot "dati_a.dat" u 1:2 w l title "Metodo 1", "dati_a.dat" u 1:3 w l title "Metodo 2"

EOFMarker