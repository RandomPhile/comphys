#!/bin/bash
g++ es2_a.cpp -o es2_a.out
./es2_a.out > dati_a.dat
pkill -x gnuplot, gnuplot_qt
#!/usr/bin/gnuplot -persist
gnuplot -persist <<-EOFMarker

set term qt font "Helvetica"
set title "Esercizio 2a"
set xlabel "h"
set ylabel "Δ_h"
set grid
set format x "%2.0tx10^{%L}"
set format y "%2.0tx10^{%L}"
set logscale x
plot "dati_a.dat" u 1:2 w l title "Differenza in avanti", \\
"dati_a.dat" u 1:3 w l title "Differenza centrale"

EOFMarker