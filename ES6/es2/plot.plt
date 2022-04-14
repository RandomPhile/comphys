set datafile sep '\t'
set grid

#set multiplot layout 1,2 
set term qt 0 title "E(t)" font "Helvetica" 
plot "dati.dat" i 0 u 1:2 w l title "E" lc "red",\
 "dati.dat" i 0 u 1:3 w l title "K" lc "blue"

set term qt 1 title "T(t)" font "Helvetica" 
plot "dati.dat" i 0 u 1:4 w l title "vel verlet" lc "red"
# set term qt 1 title "v" font "Helvetica" 
# plot "dati.dat" i 0 u 1:3 w l title "vel verlet" lc "red"

# set term qt 2 title "E" font "Helvetica" 
# plot "dati.dat" i 0 u 1:4 w l title "vel verlet" lc "red"


pause -1
