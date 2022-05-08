set datafile sep '\t'
set grid
set key left top

#set logscale xy
set format "%g"
set style line 1 lt 1 lw 1 pt 7 ps 0.7 lc rgb "red"
set style line 2 lt 1 lw 1 pt 7 ps 0.7 lc rgb "blue"
set style line 3 lt 1 lw 1 pt 7 ps 0.7 lc rgb "black"
set style line 4 lt 1 lw 1 pt 7 ps 0.7 lc rgb "orange"
set style line 5 lt 1 lw 1 pt 7 ps 0.7 lc rgb "cyan"
set style line 6 lt 1 lw 1 pt 7 ps 0.7 lc rgb "dark-grey"

if (errore == 0) {
set term qt 0 title "massa - raggio" position 0,0 font "Helvetica, 14"
set xlabel "R [km]"
set ylabel "M [masse solari]"
set key top right
set xtics add (3,10,50)
if (relativ == 1) {
	unset logscale x
	unset logscale y
} else {
	set logscale y#relativ
	set logscale x#relativ
	set xrange [3:50]#relativ
}
plot "dati.dat" i 0 u 5:6 w linespoint ls 1 title columnheader(1),\
	 "dati.dat" i 1 u 5:6 w linespoint ls 2 title columnheader(1),\
	 "dati.dat" i 2 u 5:6 w linespoint ls 3 title columnheader(1)
unset xrange
set xtics auto

if (relativ == 0) {
set term qt 1 title "raggio - pressione centrale" position 0,400 font "Helvetica, 14"
set xlabel "P_c [unità adimensionali]"
set ylabel "R [km]"
set logscale x
set logscale y
set ytics add (1,3,10,50,100)
set arrow from 4e-7,3 to 1e7,3 nohead dt "-" lc rgb "blue";
set arrow from 4e-7,50 to 1e7,50 nohead dt "-" lc rgb "blue";
plot "dati.dat" i 0 u 1:5 w linespoint ls 1 title columnheader(1),\
	 "dati.dat" i 1 u 1:5 w linespoint ls 2 title columnheader(1),\
	 "dati.dat" i 2 u 1:5 w linespoint ls 3 title columnheader(1)
set ytics auto
}
}
if (relativ == 0) {
set term qt 2 title "fattore di controllo - pressione centrale" position 640,0 font "Helvetica, 14"
set xlabel "P_c [unità adimensionali]"
set ylabel "Q [unità adimensionali]"
set key left bottom
set logscale x
unset logscale y
if (errore == 0) {
plot "dati.dat" i 0 u 1:4 w linespoint ls 1 title columnheader(1),\
	 "dati.dat" i 1 u 1:4 w linespoint ls 2 title columnheader(1),\
	 "dati.dat" i 2 u 1:4 w linespoint ls 3 title columnheader(1)
} else {
plot "dati.dat" i 0 u 1:4 w linespoint ls 1 title columnheader(1),\
	 "dati.dat" i 1 u 1:4 w linespoint ls 2 title columnheader(1),\
	 "dati.dat" i 2 u 1:4 w linespoint ls 3 title columnheader(1),\
	 "dati.dat" i 3 u 1:4 w linespoint ls 4 title "",\
	 "dati.dat" i 4 u 1:4 w linespoint ls 5 title "",\
	 "dati.dat" i 5 u 1:4 w linespoint ls 6 title ""
}
}

if (relativ == 1 && errore == 0) {
set term qt 3 title "massa - pressione centrale" position 0,400 font "Helvetica, 14"
set xlabel "P_c [unità adimensionali]"
set ylabel "M [unità adimensionali]"
set logscale x
unset logscale y
load "rect.dat"
plot "dati.dat" i 0 u 1:3 w linespoint ls 1 title columnheader(1),\
	 "dati.dat" i 1 u 1:3 w linespoint ls 2 title columnheader(1),\
	 "dati.dat" i 2 u 1:3 w linespoint ls 3 title columnheader(1)
}

pause -1