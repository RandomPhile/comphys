set datafile sep '\t'
set grid

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 0 title "blocking" font "Helvetica" position 50,0
set ytics add ("1 " 1)
#set logscale x
plot "blocking.dat" i 0 u 1:2 w lp ls 1 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 1 title "blocking2" font "Helvetica" position 50,0
set ytics add ("1 " 1)
set logscale y
set logscale x
#plot "blocking2.dat" i 0 u 1:2 w p title ""

pause-1