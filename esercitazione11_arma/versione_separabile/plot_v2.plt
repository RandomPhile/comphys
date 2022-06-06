set datafile sep '\t'
set grid

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 0 title "potenziale istant" font "Helvetica" position 400,0
set ytics add ("1 " 1)
plot "out/oss.dat" i 0 u 1:4 w l ls 1 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 1 title "Pressione" font "Helvetica" position 400,300
set ytics add ("1 " 1)
plot "out/oss.dat" i 0 u 1:5 w l ls 1 title "",\
	 		 	"" i 0 u 1:3 w l ls 1 title ""

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set term qt 4 title "errore pressione" font "Helvetica" position 50,0
set ytics add ("1 " 1)
plot "out/oss.dat" i 0 u 1:6 w l ls 1 title ""

pause-1