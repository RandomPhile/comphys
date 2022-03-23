#set term qt font "Helvetica" position 0,0
set datafile sep '\t'
set grid

#set multiplot layout 1,2
set term qt 0 title "P"
plot "dati.dat" i 0 u 1:2 w l

set term qt 1 title "m"
plot "dati.dat" i 0 u 1:3 w l


pause -1
