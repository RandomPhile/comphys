#set term qt font "Helvetica" position 0,0
set datafile sep '\t'
set grid

N = 10

set term qt 0 title "P(r)" position 0,0
set xlabel "r"
set ylabel "P"
plot "dati.dat" i 0*N u 1:2 w l lt rgb "violet" title columnheader(1),\
	 "dati.dat" i 1*N u 1:2 w l lt rgb "red" title columnheader(1), \
	 "dati.dat" i 2*N u 1:2 w l title columnheader(1)

set term qt 1 title "m(r)" position 0,0
set xlabel "r"
set ylabel "m"
plot "dati.dat" i 0*N u 1:3 w l lt rgb "violet" title columnheader(1),\
	 "dati.dat" i 1*N u 1:3 w l lt rgb "red" title columnheader(1), \
	 "dati.dat" i 2*N u 1:3 w l title columnheader(1)

set term qt 2 title "cost(Pc)" position 0,0
set xlabel "r"
set ylabel "m"
plot "dati1.dat" i 0 u 1:4 w l lt rgb "violet" title columnheader(1),\
	 "dati1.dat" i 1 u 1:4 w l lt rgb "red" title columnheader(1), \
	 "dati1.dat" i 2 u 1:4 w l title columnheader(1)

# set term qt 1 title "m" position 800,0
# plot 	"dati.dat" i 0 u 1:3 w l lt rgb "violet", \
# 	"dati.dat" i 1 u 1:3 w l lt rgb "red", \
# 	"dati.dat" i 2 u 1:3 w l

# set term qt 2 title "rho" position 400,350
# plot 	"dati.dat" i 0 u 1:4 w l lt rgb "violet", \
# 	"dati.dat" i 1 u 1:4 w l lt rgb "red", \
# 	"dati.dat" i 2 u 1:4 w l

# set term qt 3 title "M-R" position 400,350
# set logscale y
# set logscale x
# plot 	"dati.dat" i STATS_blocks-3 u 2:1  w lp lt rgb "violet", \
# 	"dati.dat" i STATS_blocks-2 u 2:1  w lp lt rgb "red", \
# 	"dati.dat" i STATS_blocks-1 u 2:1  w lp
	

#plot 'file.dat' i STATS_blocks-1 usato per accedere all'ultimo blocco

pause -1
