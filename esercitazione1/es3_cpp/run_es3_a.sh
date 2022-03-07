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
set terminal qt 0
plot "dati_a.dat" index 0 u 1:2 w l title "float", \
			   "" index 0 u 1:3 w l title "double"
set terminal qt 1
plot "dati_a.dat" index 1 u 1:2 w l title "float", \
			   "" index 1 u 1:3 w l title "double"
EOFMarker