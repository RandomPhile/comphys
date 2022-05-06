#!/bin/bash
g++ main.cpp -o main.out
time ./main.out
pkill -x gnuplot, gnuplot_qt
#gnuplot 'plot.plt'
i=0
while read line; do    
	param[$i]=$line
	i=$i+1
done < gnuplot.dat

case ${param[0]} in

	0)
		gnuplot "plot.plt"
		;;

	1)
		gnuplot -e "N=${param[1]}" -e "N_step=${param[2]}" -e "L=${param[3]}" -e "pausa=${param[4]}" -e "skip=${param[5]}" -e "dt=${param[6]}" "animation.plt"	
		;;
esac