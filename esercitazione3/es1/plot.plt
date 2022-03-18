set term qt font "Helvetica" position 0,0
set datafile sep '\t'
set grid

#set multiplot layout 1,2
set term qt 0
plot "dati.dat" i 0 u 1:2 w l
set term qt 1
plot "dati.dat" i 0 u 1:3 w l
set term qt 2
plot "dati.dat" i 0 u 1:4 w l
set term qt 3
plot "dati.dat" i 1 u 1:2 w l
set term qt 4
plot "dati.dat" i 1 u 1:3 w l
set term qt 5
plot "dati.dat" i 1 u 1:4 w l

pause -1
