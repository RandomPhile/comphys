#!/bin/bash
pkill -x gnuplot, gnuplot_qt
g++ *.cpp -o out/a.out
time ./out/a.out
