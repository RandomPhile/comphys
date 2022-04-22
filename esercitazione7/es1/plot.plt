set datafile sep '\t'
set grid

set term qt 5 title "Reticolo Iniziale" font "Helvetica" 
splot "dati.dat" i 0 u 1:2:3 w p,\

#set term qt 4 title "Reticolo dopo" font "Helvetica" 
#splot "dati.dat" i 2 u 1:2:3 w p,\

#set term qt 0 title "E(t)" font "Helvetica" 
#plot "dati.dat" i 1 u 1:2 w l title "E" lc "red",\
  #	"dati.dat" i 1 u 1:3 w l title "K" lc "blue",\
 #	"dati.dat" i 1 u 1:4 w l title "V" lc "blue"


#set term qt 1 title "T(t)" font "Helvetica" 
#plot "dati.dat" i 1 u 1:5 w l title "T" lc "red"

#set term qt 2 title "v" font "Helvetica" 
#plot "dati.dat" i 1 u 1:3 w l title "vel verlet" lc "red"

#set term qt 3 title "E" font "Helvetica" 
#plot "dati.dat" i 1 u 1:4 w l title "vel verlet" lc "red"

	 
pause-1