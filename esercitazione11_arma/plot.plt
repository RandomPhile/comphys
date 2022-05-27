set datafile sep '\t'
set grid

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 0 title "potenziale medio" font "Helvetica" position 0,0
set ytics add ("1 " 1)
plot "dati.dat" i 0 u 1:2 w l ls 1 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 1 title "Pressione scalata media" font "Helvetica" position 0,400
set ytics add ("1 " 1)
plot "dati.dat" i 0 u 1:3 w l ls 1 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 2 title "potenziale istant" font "Helvetica" position 400,0
set ytics add ("1 " 1)
plot "dati.dat" i 0 u 1:4 w l ls 1 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 3 title "Pressione scalata istant" font "Helvetica" position 400,400
set ytics add ("1 " 1)
plot "dati.dat" i 0 u 1:5 w l ls 1 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 4 title "g(r)" font "Helvetica" position 50,0
set ytics add ("1 " 1)
plot "dati.dat" i 1 u 1:2 w l ls 1 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 5 title "reticolo" font "Helvetica" position 800,0
set ytics add ("1 " 1)
splot "dati.dat" i 2 u 1:2:3 w p 

# set term qt 0 title "Hist" font "Helvetica" position 0,350
# plot "dati.dat" i 0 u 1:2 w l title "freq" lc "red"

pause-1