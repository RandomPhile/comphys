set term qt font "Helvetica" position 0,0
set datafile sep '\t'
set grid

set term qt 0 title "Reticolo Iniziale" font "Helvetica" position 0,0
plot "dati.dat" i 0 u 1:2 w p

pause -1
