#!/bin/bash
g++ main.cpp -o main.out
./main.out
pkill -x gnuplot, gnuplot_qt
gnuplot 'plot.plt'