set datafile sep '\t'
set grid

set term qt 0 title "g(r)" font "Helvetica" position 0,350
plot "dati.dat" i 1 u 1:2 w lp title "" lc "red"

set term qt 1 title "reticolo" font "Helvetica" position 0,0
splot "dati.dat" i 0 u 1:2:3 w p title "" lc "red"

pause-1