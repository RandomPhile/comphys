set datafile sep '\t'
set grid

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5

set term qt 0 title "Pressione scalata - densità" font "Helvetica" position 0,0
set xlabel "{/Symbol r}" font ",13"
set ylabel "Q" font ",13"
set ytics add ("1 " 1)
plot "out/dati.dat" i 0 u 1:2:3 notitle w yerrorlines lw 1.5 lc rgb '#0060ad' 

pause-1