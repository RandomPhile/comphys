#!/bin/bash
g++ main.cpp -o main -std=c++11 -O2 -larmadillo -o main.out
time ./main.out
pkill -x gnuplot, gnuplot_qt
gnuplot 'plot.plt'