set datafile sep '\t'
set grid

#set multiplot layout 1,2 
set term qt 0 title "x" font "Helvetica" 
plot "dati.dat" i 0 u 1:2 w l title "eulero exp",\
	 "dati.dat" i 1 u 1:2 w l title "eulero cromer",\
	 "dati.dat" i 2 u 1:2 w l title "eulero imp",\
	 "dati.dat" i 3 u 1:2 w l title "runge kutta",\
	 "dati.dat" i 4 u 1:2 w l title "vel verlet" lc "red"

set term qt 1 title "v" font "Helvetica" 
plot "dati.dat" i 0 u 1:3 w l title "eulero exp",\
	 "dati.dat" i 1 u 1:3 w l title "eulero cromer",\
	 "dati.dat" i 2 u 1:3 w l title "eulero imp",\
	 "dati.dat" i 3 u 1:3 w l title "runge kutta",\
	 "dati.dat" i 4 u 1:3 w l title "vel verlet" lc "red"

set term qt 2 title "E" font "Helvetica" 
plot "dati.dat" i 0 u 1:4 w l title "eulero exp",\
	 "dati.dat" i 1 u 1:4 w l title "eulero cromer",\
	 "dati.dat" i 2 u 1:4 w l title "eulero imp",\
	 "dati.dat" i 3 u 1:4 w l title "runge kutta",\
	 "dati.dat" i 4 u 1:4 w l title "vel verlet" lc "red"


pause -1
