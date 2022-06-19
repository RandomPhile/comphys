set datafile sep '\t'
set grid

set term qt 0 title "K, V, E" font "Helvetica" position 0,0
set xlabel "Passi temporali"
set xlabel font ",13"
set key center left
set key font ",12"
plot "out/osservabili.dat" i 0 u 1:2 w l title "K/{/Symbol e}" lc "red",\
  	 "out/osservabili.dat" i 0 u 1:3 w l title "V/{/Symbol e}" lc "blue",\
 	 "out/osservabili.dat" i 0 u 1:4 w l title "E/{/Symbol e}" lc "green"
   
   
set term qt 1 title "T(t)" font "Helvetica" position 640,0
set xlabel "Passi temporali"
set xlabel font ",13"
set ylabel "T/T_0"
set ylabel font ",13"
plot "out/osservabili.dat" i 0 u 1:5 w l title "" lc "red"



set term qt 2 title "Q(t)" font "Helvetica" position 0,350
set xlabel "Passi temporali"
set xlabel font ",13"
set ylabel "Q"
set ylabel font ",13"

plot "out/osservabili.dat" i 0 u 1:6 w l title "" lc "red"

# set term qt 0 title "Hist" font "Helvetica" position 0,350
# plot "dati.dat" i 0 u 1:2 w l title "freq" lc "red"

pause-1
