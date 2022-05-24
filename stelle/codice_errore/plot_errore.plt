set logscale x
set xlabel "h"
set ylabel "{/Symbol e}"
set xlabel font ",14"
set ylabel font ",14"
set grid

set grid


plot "err.txt" w p pointtype 7 linecolor rgb "red" title "Errore relativo (Eulero)" #EULERO
#plot "err.txt" w p pointtype 7 linecolor rgb "black" title "Errore relativo (Runge-Kutta)" #RK



pause -1


