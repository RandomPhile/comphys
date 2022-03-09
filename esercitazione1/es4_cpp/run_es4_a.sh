#!/bin/bash
g++ es4_a.cpp -o es4_a.out
./es4_a.out
pkill -x gnuplot, gnuplot_qt
#!/usr/bin/gnuplot -persist
gnuplot -persist <<-EOFMarker

set term qt font "Helvetica" #position 0,0
#set term eps font "Helvetica"
#set output "es4a.eps"

set datafile sep '\t'
#set title "Esercizio 4a"
set xlabel "x0"
set ylabel "ε_r"
set grid
set logscale x
#set xrange [1e-20:1e-14]
set format x "%2.0tx10^{%L}"

set multiplot layout 2,1
#set xrange [1e-20:1e-1]
plot "dati1.dat" u 1:2 w l title "x_i = 1.0",\\
"dati2.dat" u 1:2 w l title "x_i = 0.1"

#set yrange [0:6e-8]
plot "dati1.dat" u 1:3 w l title "x_i = 1.0 new",\\
"dati2.dat" u 1:3 w l title "x_i = 0.1 new"

#, "dati1.dat" u 1:3 w l title "2"

#plot "dati2.dat" u 1:2 w l title "x_i = 0.1"
#, "dati2.dat" u 1:3 w l title "2"

EOFMarker