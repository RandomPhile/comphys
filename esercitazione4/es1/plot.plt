#set term qt font "Helvetica" position 0,0
set datafile sep '\t'
set grid

#set multiplot layout 1,2
set term qt 0 title "P" position 0,0
plot 	"dati.dat" i 0 u 1:2 w l, \
	"dati.dat" i 1 u 1:2 w l, \
	"dati.dat" i 2 u 1:2 w l

set term qt 1 title "m" position 800,0
plot 	"dati.dat" i 0 u 1:3 w l, \
	"dati.dat" i 1 u 1:3 w l, \
	"dati.dat" i 2 u 1:3 w l

set term qt 2 title "rho" position 400,350
plot 	"dati.dat" i 0 u 1:4 w l, \
	"dati.dat" i 1 u 1:4 w l, \
	"dati.dat" i 2 u 1:4 w l

pause -1
