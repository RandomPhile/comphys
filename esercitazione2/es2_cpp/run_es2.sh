#!/bin/bash
g++ es2.cpp -o es2.out
./es2.out
pkill -x gnuplot, gnuplot_qt

#!/usr/bin/gnuplot -persist
gnuplot -persist <<-EOFMarker

set term qt font "Helvetica" #position 0,0
set datafile sep '\t'
set grid

f(e,v) = 1/tan(sqrt(v-e))
g(e,v) = -sqrt(e/(v-e))
#set xrange[0.15:0.2]
#plot f(x,3.45804), g(x,3.45804)

set multiplot layout 1,2
plot "dati1.dat" u 1:2 w l
plot "dati2.dat" u 1:2 w l
EOFMarker