set datafile sep '\t'
set grid

set term qt 0 title "Gaussian?" font "Helvetica" 
plot "dati.dat" i 0 u 2:1 w p,\

	 
pause-1