#!/bin/bash
g++ es2_b.cpp -o es2_b.out
./es2_b.out > dati_b.dat
pkill -x gnuplot, gnuplot_qt
#!/usr/bin/gnuplot -persist
gnuplot -persist <<-EOFMarker

set term qt font "Helvetica"
set title "Esercizio 2b"
set xlabel "h"
set ylabel "Δ_h"
set grid
set format x "%2.0tx10^{%L}"
set format y "%2.0tx10^{%L}"
set logscale y
set logscale x

plot "dati_b.dat" u 1:2 w l title "Differenza in avanti", \\
"dati_b.dat" u 1:3 w l title "Differenza centrale (senza errore su 2*h)", \\
"dati_b.dat" u 1:4 w l title "Differenza centrale (con errore su 2*h)"

EOFMarker
