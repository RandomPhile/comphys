set datafile sep '\t'
set grid

set term qt 0 title "g(r)" font "Helvetica" position 0,350
plot "dati.dat" i 0 u 1:2 w lp title "" lc "red"

pause-1