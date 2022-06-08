set datafile sep '\t'
set grid

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 1 pi -1 ps 1.5
set term qt 0 title "blocking Pressione" font "Helvetica" position 50,0
set ytics add ("1 " 1)
plot "blocking.dat" i 0 u 1:2 w p ls 1 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 1 pi -1 ps 1.5
set term qt 1 title "blocking Energia" font "Helvetica" position 450,0
set ytics add ("1 " 1)
plot "blocking.dat" i 0 u 1:3 w p ls 1 title ""

pause-1