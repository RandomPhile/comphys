set datafile sep '\t'
set grid
# set key center left
set format "%g"
set style line 1 lc rgb "red"
set style line 2 lc rgb "blue"
set style line 3 lc rgb "green"

set term qt 0 title "Energia" position 0,0 font "Helvetica, 14"
set xlabel "passi temporali"
plot "dati.dat" i 0 u 1:2 w l ls 1 title "K/ε",\
	 "dati.dat" i 0 u 1:3 w l ls 2 title "V/ε",\
	 "dati.dat" i 0 u 1:4 w l ls 3 title "E/ε"

set term qt 1 title "Pressione" position 0,300 font "Helvetica, 14"
set xlabel "Passi temporali"
plot "dati.dat" i 0 u 1:5 w l ls 1 title "Pressione istantanea",\
	 "dati.dat" i 0 u 1:6 w l ls 2 title "Pressione media"

pause -1