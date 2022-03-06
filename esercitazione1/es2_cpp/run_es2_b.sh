#!/bin/bash
g++ es2_b.cpp -o es2_b.out
./es2_b.out > dati_b.dat
pkill -x gnuplot, gnuplot_qt
#!/usr/bin/gnuplot -persist
gnuplot -persist <<-EOFMarker

set termoption enhanced
set title "Esercizio 2b"
set xlabel "h"
set ylabel "\{/Symbol D}_h"
set grid
set logscale x
plot "dati_b.dat" u 1:2 w l title "Metodo 1", "dati_b.dat" u 1:3 w l title "Metodo 2", "dati_b.dat" u 1:4 w l title "Metodo 3"

EOFMarker