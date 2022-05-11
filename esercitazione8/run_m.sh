export CPATH=/opt/homebrew/include
export LIBRARY_PATH=/opt/homebrew/lib

#!/bin/bash
g++ main_arma.cpp -o main_arma -std=c++11 -O2 -larmadillo -o main.out
./main.out
pkill -x gnuplot, gnuplot_qt
gnuplot 'plot.plt'