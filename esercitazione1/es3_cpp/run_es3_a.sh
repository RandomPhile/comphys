#!/bin/bash
g++ es3_a.cpp -o es3_a.out
./es3_a.out
pkill -x gnuplot, gnuplot_qt
#!/usr/bin/gnuplot -persist
gnuplot -persist <<-EOFMarker

set termoption enhanced
set datafile sep ','
set term qt font "Arial"
#set title "Esercizio 3a"
#set xlabel "h"
#set ylabel "\{/Symbol D}_h"
set grid
set logscale x

set terminal qt 0
plot "dati1.dat" u 1:2 w l title "x1"
set terminal qt 1
plot "dati2.dat" u 1:2 w l title "x1s"
set terminal qt 2
plot "dati1.dat" u 1:3 w l title "x2"
set terminal qt 3
plot "dati2.dat" u 1:3 w l title "x2s"

EOFMarker