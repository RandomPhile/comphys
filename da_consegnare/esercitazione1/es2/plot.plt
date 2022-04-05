set term qt font "Helvetica,12" position 0,0

set xlabel "h"
set ylabel "Δ_h"
set grid
set format x "%2.0tx10^{%L}"
set format y "%2.0tx10^{%L}"
set logscale x
set logscale y
set multiplot layout 2,1
set title "Differenza in avanti"
plot "dati.dat" i 0 u 1:2 w l title "double",\
			 "" i 1 u 1:2 w l title "float" lt rgb "red"
set title "Differenza centrale"
plot "" i 0 u 1:3 w l title "double",\
     "" i 1 u 1:3 w l title "float" lt rgb "red" 

pause -1
