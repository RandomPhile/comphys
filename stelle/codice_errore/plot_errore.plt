set logscale x
set xlabel "h"
set ylabel "{/Symbol e}"
set xlabel font ",14"
set ylabel font ",14"
set grid

set grid

#set xrange [9*1e-8:6*1e-4] #E
#set yrange [-0.005:0.1] #E

set xrange [1e-5:6*1e-4] #RK
set yrange [-0.5*1e-5:1e-4] #RK

set key top left
set key font ",11"
#plot "err.txt" w p pointtype 7 linecolor rgb "red" title "Errore relativo (Eulero)" #E
plot "err.txt" w p pointtype 7 linecolor rgb "black" title "Errore relativo (Runge-Kutta)" #RK



pause -1


