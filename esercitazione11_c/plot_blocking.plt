set datafile sep '\t'
set grid
# set key center left
set format "%g"
set style line 1 lt 1 lw 1 pt 7 ps 0.7 lc rgb "red"
set style line 2 lt 1 lw 1 pt 7 ps 0.7 lc rgb "blue"
set style line 3 lt 1 lw 1 pt 7 ps 0.7 lc rgb "black"

set term qt 0 title "ΔP, ΔE" position 0,0 font "Helvetica, 14"
set xlabel "B"
# set logscale x
plot "dati_blocking.dat" i 0 u 1:2 w linespoint ls 1 title "ΔP",\
	 "dati_blocking.dat" i 0 u 1:3 w linespoint ls 1 title "ΔE"

pause -1