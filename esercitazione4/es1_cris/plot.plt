set term qt 0

set xlabel "r (km)"
set ylabel "P (MeV/fm^3)"
set xlabel font ",12"
set ylabel font ",12"
set grid
set key top right
set key font ",10"
plot "stella10.txt" using 1:2 with lines linecolor rgb "red" title "Modello con {/Symbol G} = 5/3"


set term qt 1

set xlabel "r (km)"
set ylabel "m (unità di masse solari)"
set xlabel font ",12"
set ylabel font ",12"
set grid
set key top left
set key font ",10"
plot "stella10.txt" using 1:3 with lines linecolor rgb "red" title "Modello con {/Symbol G} = 5/3"

set term qt 2

set xlabel "R (km)"
set ylabel "M (unità di masse solari)"
set xlabel font ",12"
set ylabel font ",12"
set grid
set key top right
set key font ",10"
plot "massa_raggio.txt" using 1:2 index 0 with lines linecolor rgb "red" title "{/Symbol G} = 5/3" #points pointtype 7 
replot "massa_raggio.txt" using 1:2 index 1 with lines linecolor rgb "blue" title "{/Symbol G} = 4/3"
replot "massa_raggio.txt" using 1:2 index 2 with lines linecolor rgb "black" title "{/Symbol G} = 2.54"

set term qt 3

set xlabel "R (km)"
set ylabel "Fattore di controllo (unità adimensionali)"
set xlabel font ",12"
set ylabel font ",12"
set grid
set key top right
set key font ",10"
plot "massa_raggio.txt" using 1:3 index 0 with lines linecolor rgb "red" title "Fattore per {/Symbol G} = 5/3" #points pointtype 7 
replot "massa_raggio.txt" using 1:3 index 1 with lines linecolor rgb "blue" title "Fattore per {/Symbol G} = 4/3"
replot "massa_raggio.txt" using 1:3 index 2 with lines linecolor rgb "black" title "Fattore per {/Symbol G} = 2.54"

pause -1


