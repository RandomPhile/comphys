set logscale x
set xlabel "h"
set ylabel "{/Symbol e}"
set xlabel font ",14"
set ylabel font ",14"
set grid
set key at 3*1e-5, 0.142 #EULERO
#set key at 0.00023, 1.96*1e-5 #RK
set key font ",11"
set grid
set yrange [-0.005:0.15]   #EULERO
set xrange [9*1e-7:1.5*1e-3]   #EULERO
#set yrange [-0.000001:0.00002] #RK
#set xrange [1.5*1e-5: 2*1e-3]  #RK

plot "err.txt" w p pointtype 7 linecolor rgb "red" title "Errore relativo (Eulero)" #EULERO
#plot "err.txt" w p pointtype 7 linecolor rgb "black" title "Errore relativo (Runge-Kutta)" #RK



pause -1


