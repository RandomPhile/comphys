set term qt font "Helvetica" position 0,0

set title "Esercizio 2b"
set xlabel "h"
set ylabel "Δ_h"
set grid
set format x "%2.0tx10^{%L}"
set format y "%2.0tx10^{%L}"
set logscale y
set logscale x

plot "dati.dat" u 1:2 w l title "Differenza in avanti", \
			 "" u 1:3 w l title "Differenza centrale (senza errore su 2*h)", \
			 "" u 1:4 w l title "Differenza centrale (con errore su 2*h)"

pause -1
