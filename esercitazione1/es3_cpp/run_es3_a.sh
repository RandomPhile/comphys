#!/bin/bash
g++ es3_a.cpp -o es3_a.out
./es3_a.out > dati_a.dat
pkill -x gnuplot, gnuplot_qt
#!/usr/bin/gnuplot -persist
gnuplot -persist <<-EOFMarker

set termoption enhanced
set title "Esercizio 3a"
set xlabel "h"
set ylabel "\{/Symbol D}_h"
set grid
set logscale x
plot "dati_a.dat" u 1:2 w l title "Metodo 1", "dati_a.dat" u 1:3 w l title "Metodo 2"

EOFMarker