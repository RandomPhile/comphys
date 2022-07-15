set datafile sep '\t'
set grid

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 1 pi -1 ps 1.5
set term qt 0 title "Errore Pressione - #passi per blocco blocking" font "Helvetica" position 50,0
set ytics add ("1 " 1)
set xlabel "B" font ",14"
set ylabel "{/Symbol s}_{P}" font ",14"
plot "out/blocking.dat" i 0 u 1:2 w lp ls 1 title ""


set term qt 1 title "Errore Energia - #passi per blocco blocking" font "Helvetica" position 750,0
set ytics add ("1 " 1)
set xlabel "B" font ",14"
set ylabel "{/Symbol s}_{E}" font ",14"
plot "out/blocking.dat" i 0 u 1:3 w lp ls 1 title ""


set term qt 2 title "Errore Pressione - #passi per blocco jackknife" font "Helvetica" position 50,350
set ytics add ("1 " 1)
set xlabel "B" font ",14"
set ylabel "{/Symbol s}_{P}" font ",14"
plot "out/jackknife.dat" i 0 u 1:2 w lp ls 1 title ""


set term qt 3 title "Errore Energia - #passi per blocco jackknife" font "Helvetica" position 750,350
set ytics add ("1 " 1)
set xlabel "B" font ",14"
set ylabel "{/Symbol s}_{E}" font ",14"
plot "out/jackknife.dat" i 0 u 1:3 w lp ls 1 title ""

pause-1