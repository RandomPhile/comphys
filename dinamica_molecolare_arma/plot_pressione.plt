set datafile sep '\t'
set grid

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5

set term qt 0 title "Pressione scalata - densità" font "Helvetica" position 0,0
set xlabel "{/Symbol r}"
set ylabel "P / ({/Symbol r} k_B T_{req})"
set ytics add ("1 " 1)
plot "out/pressione.dat" i 0 u 1:2 w linespoints ls 1 title ""


# set term qt 0 title "Hist" font "Helvetica" position 0,350
# plot "dati.dat" i 0 u 1:2 w l title "freq" lc "red"

pause-1