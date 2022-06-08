#!/bin/bash
pkill -x gnuplot, gnuplot_qt
gcc main.c -o main.out
time ./main.out
