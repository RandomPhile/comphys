set datafile sep '\t'
set grid

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 0 title "P/T_req (t)" font "Helvetica" position 0,0
plot "dati.dat" i 0 u 1:2 w lp title "" lc "blue"

pause-1