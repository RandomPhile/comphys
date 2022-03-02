#!/bin/bash
gcc es2.c -o es2
./es2 > dati
pkill -x gnuplot, gnuplot_qt
#!/usr/bin/gnuplot -persist
gnuplot -persist <<-EOFMarker
	set logscale x
	plot "dati" u 1:2
EOFMarker