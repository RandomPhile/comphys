set term qt font "Helvetica,15" position 0,0
set datafile sep '\t'
set grid
set xlabel "e"
set ylabel "f(e)"
set xrange [0:1]
f(e,v) = 1/tan(sqrt(v-e))
g(e,v) = -sqrt(e/(v-e))
plot f(x,3.45804) title "LHS", g(x,3.45804) title "RHS"
#plot f(x,2.461047) title "LHS", g(x,2.461047) title "RHS"

pause -1
