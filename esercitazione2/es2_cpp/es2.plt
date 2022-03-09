set term qt font "Helvetica" position 0,0
set datafile sep '\t'
set grid

f(e,v) = 1/tan(sqrt(v-e))
g(e,v) = -sqrt(e/(v-e))
#plot f(x,3.45804), g(x,3.45804)
#plot f(x,2.461047), g(x,2.461047)

#set multiplot layout 1,2
plot "dati.dat" i 0 u 1:2 w l title columnheader(1), \
			 "" i 1 u 1:2 w l title columnheader(1), \
			 "" i 2 u 1:2 w l title columnheader(1)
pause -1
