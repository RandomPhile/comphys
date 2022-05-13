#!/bin/bash
g++ main.cpp -o main.out
time ./main.out
pkill -x gnuplot, gnuplot_qt

gnuplot "plot.plt"

# i=0
# while read line; do    
# 	param[$i]=$line
# 	i=$i+1
# done < gnuplot.dat

# if ((${param[0]} == 0));
# then
# 	gnuplot "plot.plt"
# else
# 	gnuplot "plot2.plt"
# fi