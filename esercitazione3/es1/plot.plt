set term qt font "Helvetica" position 0,0
set datafile sep '\t'
set grid


set term qt 0
set title "x Eulero esplicito"
plot "dati.dat" i 0 u 1:2 w l
set term qt 1
set title "v Eulero esplicito"
plot "dati.dat" i 0 u 1:3 w l
set term qt 2
set title "E Eulero esplicito"
plot "dati.dat" i 0 u 1:4 w l

set term qt 3
set title "x Eulero Cromer"
plot "dati.dat" i 1 u 1:2 w l
set term qt 4
set title "v Eulero Cromer"
plot "dati.dat" i 1 u 1:3 w l
set term qt 5
set title "E Eulero Cromer"
plot "dati.dat" i 1 u 1:4 w l

set term qt 6
set title "x Eulero implicito"
plot "dati.dat" i 2 u 1:2 w l
set term qt 7
set title "v Eulero implicito"
plot "dati.dat" i 2 u 1:3 w l
set term qt 8
set title "E Eulero implicito"
plot "dati.dat" i 2 u 1:4 w l

pause -1
