set term qt font "Helvetica" position 0,0

set datafile sep '\t'
set xlabel "b"
set ylabel "ε_r"
set grid
set logscale x

set multiplot layout 2,2
plot "dati1.dat" u 1:2 w l title "x1"
plot "dati1.dat" u 1:3 w l title "x2"
plot "dati2.dat" u 1:2 w l title "x1s"
plot "dati2.dat" u 1:3 w l title "x2s"

pause -1
