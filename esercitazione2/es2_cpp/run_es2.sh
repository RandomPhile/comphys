#!/bin/bash
g++ es2.cpp -o es2.out
./es2.out
pkill -x gnuplot, gnuplot_qt
gnuplot 'es2.plt'