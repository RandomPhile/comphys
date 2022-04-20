set datafile sep '\t'
set grid

set term qt 0 title "Reticolo Iniziale" font "Helvetica" 
splot "dati.dat" i 0 u 1:2:3 w p,\

	 
pause-1