set term qt font "Helvetica" position 0,0
#set term eps font "Helvetica"
#set output "es4a.eps"
set datafile sep '\t'
set grid

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
set yrange [0:1e-7]
plot "dati1.dat" u 1:2 w l title "x_i = 1.0",\
	 "dati2.dat" u 1:2 w l title "x_i = 0.1"

set yrange [0:1e-7]
#set yrange [0:6e-8]
plot "dati1.dat" u 1:3 w l title "x_i = 1.0 new",\
	 "dati2.dat" u 1:3 w l title "x_i = 0.1 new"

pause -1
