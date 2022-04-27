set datafile sep '\t'
set grid
set key left top

#set logscale xy
set format "%g"

set term qt 0 title "massa - raggio" position 0,0
#set xlabel "R [unità adimensionali]"#u 2:3
set xlabel "R x R_0    [km]"
set ylabel "M(R) [unità adimensionali]"
#set xrange [0:3]
unset logscale y
unset logscale x
plot "dati.dat" i 0 u 1:3 w l lt rgb "red"   title columnheader(1),\
	 "dati.dat" i 1 u 1:3 w l lt rgb "blue"  title columnheader(1), \
	 "dati.dat" i 2 u 1:3 w l lt rgb "black" title columnheader(1)
unset xrange
set logscale x

# set term qt 1 title "fattore di controllo - pressione centrale" position 0,0
# set xlabel "P_c [unità adimensionali]"
# set ylabel "cost(P_c) [unità adimensionali]"
# unset logscale y
# plot "dati.dat" i 0 u 5:4 w l lt rgb "red"   title columnheader(1),\
# 	 "dati.dat" i 1 u 5:4 w l lt rgb "blue"  title columnheader(1), \
# 	 "dati.dat" i 2 u 5:4 w l lt rgb "black" title columnheader(1)
# set logscale y

# set term qt 3 title "raggio - pressione centrale" position 0,0
# set xlabel "P_c"
# set ylabel "R(P_c) x R_0    [km]"
# set ytics (1,3,10,50,100)
# set arrow from 4e-7,3 to 1e7,3 nohead dt "-" lc rgb "blue";
# set arrow from 4e-7,50 to 1e7,50 nohead dt "-" lc rgb "blue";
# plot "dati.dat" i 0 u 1:5 w l lt rgb "red"   title columnheader(1),\
# 	 "dati.dat" i 1 u 1:5 w l lt rgb "blue"  title columnheader(1), \
# 	 "dati.dat" i 2 u 1:5 w l lt rgb "black" title columnheader(1)

pause -1