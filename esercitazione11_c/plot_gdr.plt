set datafile sep '\t'
set grid
# set key center left
set format "%g"
set style line 1 lt 1 lw 1 pt 7 ps 0.7 lc rgb "red"
set style line 2 lt 1 lw 1 pt 7 ps 0.7 lc rgb "blue"
set style line 3 lt 1 lw 1 pt 7 ps 0.7 lc rgb "black"

set term qt 0 title "g(r)" position 0,0 font "Helvetica, 14"
set xlabel "r_k"
plot "dati_gdr.dat" i 0 u 1:3 w linespoint ls 1 title "g"#,\
	# "dati_gdr.dat" i 0 u 1:2 w linespoint ls 2 title "freq"

pause -1