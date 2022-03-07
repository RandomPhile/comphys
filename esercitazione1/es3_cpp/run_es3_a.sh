#!/bin/bash
g++ es3_a.cpp -o es3_a.out
./es3_a.out
pkill -x gnuplot, gnuplot_qt
#!/usr/bin/gnuplot -persist
gnuplot -persist <<-EOFMarker

set term qt font "Helvetica"
#set term eps font "Helvetica"
###
#set output "b_pos.eps"
#set xrange [1e2:1e5]
#set xtics ("10^2" 1e2,"10^3" 1e3,"10^4" 1e4,"10^5" 1e5)
###
#set output "b_neg.eps"
#set xrange [1e5:1e2]
#set xtics ("-10^2" 1e2,"-10^3" 1e3,"-10^4" 1e4,"-10^5" 1e5)


set datafile sep ','
#set title "Esercizio 3a"
set xlabel "b"
set ylabel "ε_r"
set grid
set logscale x


set multiplot layout 2,2
plot "dati1.dat" u 1:2 w l title "x1"
plot "dati1.dat" u 1:3 w l title "x2"
plot "dati2.dat" u 1:2 w l title "x1s"
plot "dati2.dat" u 1:3 w l title "x2s"


EOFMarker