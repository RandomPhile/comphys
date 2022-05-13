set datafile sep '\t'
set grid

# set term qt 0 title "K, V, E" font "Helvetica" position 640,0
# #set yrange [0:]
# plot "dati.dat" i 0 u 1:2 w l title "K" lc "red",\
#   	 "dati.dat" i 0 u 1:3 w l title "V" lc "blue",\
#  	 "dati.dat" i 0 u 1:4 w l title "E" lc "green"
# set yrange []

set term qt 1 title "T(t)" font "Helvetica" position 0,350
plot "dati.dat" i 0 u 1:5 w l title "" lc "red"

set term qt 2 title "P/T_req (t)" font "Helvetica" position 0,0
plot "dati.dat" i 0 u 1:6 w l title "" lc "red"

pause-1