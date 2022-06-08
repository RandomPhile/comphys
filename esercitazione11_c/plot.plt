set font "Times"
set term qt 0# title "K, V, E" font "Helvetica" position 0,0
set datafile sep '\t'
set grid
set xlabel "Passi temporali"
set xlabel font ",13"
# set key center left
set key font ",12"
plot "dati.dat" i 0 u 1:2 w l title "K/{/Symbol e}" lc "red",\
	 "dati.dat" i 0 u 1:3 w l title "V/{/Symbol e}" lc "blue",\
	 "dati.dat" i 0 u 1:4 w l title "E/{/Symbol e}" lc "green"


set term qt 1# title "Pressione" font "Helvetica" position 0,0
set xlabel "Passi temporali"
set xlabel font ",13"
# set key center left
set key font ",12"
plot "dati.dat" i 0 u 1:5 w l title "Pressione istantanea" lc "red",\
		 "dati.dat" i 0 u 1:6 w l title "Pressione media" lc "blue"

pause -1