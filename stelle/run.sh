#!/bin/bash
g++ main.cpp -o main.out
./main.out
pkill -x gnuplot, gnuplot_qt
#gnuplot 'plot.plt'
i=0
while read line; do    
	param[$i]=$line
	i=$i+1
done < gnuplot.dat

gnuplot -e "relativ=${param[0]}" -e "errore=${param[1]}" "plot.plt"