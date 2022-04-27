set datafile sep '\t'
set grid

set term qt 0 title "Reticolo Iniziale" font "Helvetica" position 0,0
splot "dati.dat" i 0 u 1:2:3 w p

set term qt 1 title "E(t)" font "Helvetica" position 800,0
plot "dati.dat" i 1 u 1:2 w l title "E" lc "red",\
  	"dati.dat" i 1 u 1:3 w l title "K" lc "blue",\
 	"dati.dat" i 1 u 1:4 w l title "V" lc "violet"


set term qt 2 title "T(t)" font "Helvetica" position 800,350
plot "dati.dat" i 1 u 1:5 w l title "T" lc "red"

set term qt 3 title "V" font "Helvetica" position 400,175
plot "dati.dat" i 1 u 1:4 w l title "V" lc "red"

set term qt 4 title "Reticolo dopo" font "Helvetica" position 0,350
splot "dati.dat" i 2 u 1:2:3 w p
	 
pause-1