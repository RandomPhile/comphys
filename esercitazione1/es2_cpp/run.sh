#!/bin/bash
g++ es2.cpp -o es2.out
./es2.out > dati.dat
pkill -x gnuplot, gnuplot_qt
#!/usr/bin/gnuplot -persist
gnuplot -persist <<-EOFMarker
	set logscale x
	plot "dati.dat" u 1:2 w linespoints
EOFMarker