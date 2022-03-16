set term qt font "Helvetica" position 0,0
set datafile sep '\t'
set grid

#set multiplot layout 1,2
plot "dati.dat" u 1:2 w l
#plot "dati.dat" u 1:3 w l
pause -1
