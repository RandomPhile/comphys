set term qt font "Helvetica,12" position 0,0
set datafile sep '\t'
set grid
set format x "%2.0tx10^{%L}"
N = 20

set term qt 0 title "P(r)" position 0,0
set xlabel "r"
set ylabel "P(r)"
set logscale y
plot "dati.dat" i 0*N u 1:2 w l lt rgb "violet" title columnheader(1),\
	 "dati.dat" i 1*N u 1:2 w l lt rgb "red" title columnheader(1), \
	 "dati.dat" i 2*N u 1:2 w l title columnheader(1)

set term qt 1 title "m(r)" position 0,0
set xlabel "r"
set ylabel "m(r)"
set logscale x
set logscale y
plot "dati.dat" i 0*N u 1:3 w l lt rgb "violet" title columnheader(1),\
	 "dati.dat" i 1*N u 1:3 w l lt rgb "red" title columnheader(1), \
	 "dati.dat" i 2*N u 1:3 w l title columnheader(1)

set term qt 2 title "cost(P_c)" position 0,0
set xlabel "P_c"
set ylabel "cost(P_c)"
set logscale x
plot "dati1.dat" i 0 u 1:4 w l lt rgb "violet" title columnheader(1),\
	 "dati1.dat" i 1 u 1:4 w l lt rgb "red" title columnheader(1), \
	 "dati1.dat" i 2 u 1:4 w l title columnheader(1)

set term qt 3 title "m(P_c)" position 0,0
set xlabel "P_c"
set ylabel "m(P_c)"
set logscale x
set logscale y
plot "dati1.dat" i 0 u 1:3 w l lt rgb "violet" title columnheader(1),\
	 "dati1.dat" i 1 u 1:3 w l lt rgb "red" title columnheader(1), \
	 "dati1.dat" i 2 u 1:3 w l title columnheader(1)

set term qt 4 title "r(P_c)" position 0,0
set xlabel "P_c"
set ylabel "r(P_c) x R_0    [km]"
set logscale x
set ytics (1,3,10,50,100)
set arrow from 4e-7,3 to 1e7,3 nohead dt "-" lc rgb "blue";
set arrow from 4e-7,50 to 1e7,50 nohead dt "-" lc rgb "blue";
plot "dati1.dat" i 0 u 1:2 w l lt rgb "violet" title columnheader(1),\
	 "dati1.dat" i 1 u 1:2 w l lt rgb "red" title columnheader(1),\
	 "dati1.dat" i 2 u 1:2 w l title columnheader(1)

pause -1
