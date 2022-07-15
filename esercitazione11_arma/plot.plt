set datafile sep '\t'
set grid

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 0 title "energia" font "Helvetica" position 400,0
set ytics add ("1 " 1)
set xlabel "Numero di passi" font ",14"
set ylabel "E" font ",14"
plot "out/dati.dat" i 0 u 1:8 w l ls 1 title "energia istantanea",\
	 		 "" i 0 u 1:4 w l lc "red" title "energia media"

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 1 title "Pressione" font "Helvetica" position 400,300
set ytics add ("1 " 1)
set xlabel "Numero di passi" font ",14"
set ylabel "P/({/Symbol r}T)" font ",14"
plot "out/dati.dat" i 0 u 1:6 w l ls 1 title "pressione istantanea",\
	 		 "" i 0 u 1:3 w l lc "red" title "pressione media"

#set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
#set term qt 3 title "errore pressione" font "Helvetica" position 50,0
#set ytics add ("1 " 1)
#plot "out/dati.dat" i 0 u 1:7 w l ls 1 title ""
#
#set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
#set term qt 4 title "errore energia" font "Helvetica" position 50,0
#set ytics add ("1 " 1)
#plot "out/dati.dat" i 0 u 1:9 w l ls 1 title ""

pause-1 