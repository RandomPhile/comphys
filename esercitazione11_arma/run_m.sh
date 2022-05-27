export CPATH=/opt/homebrew/include
export LIBRARY_PATH=/opt/homebrew/lib

#!/bin/bash
g++ main.cpp -o main -std=c++11 -O2 -larmadillo -o main.out
time ./main.out
pkill -x gnuplot, gnuplot_qt

i=0
while read line; do    
	param[$i]=$line
	i=$i+1
done < gnuplot.dat

if ((${param[0]} == -1));
then
	gnuplot "plot2.plt"
else
	gnuplot "plot.plt"
fi