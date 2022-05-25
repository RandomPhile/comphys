#!/bin/bash
pkill -x gnuplot, gnuplot_qt
g++ main.cpp -o out/main.out
time out/main.out