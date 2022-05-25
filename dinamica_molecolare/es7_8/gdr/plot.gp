set datafile sep '\t'
set grid

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5

set term qt 0 title "g(r)" font "Helvetica" position 0,0
set xlabel "r"
set ylabel "g"
# set ytics add ("1 " 1)
plot "dati.dat" i 0 u 1:3 w linespoints ls 1 title ""

pause-1