#!/bin/bash

gcc stelle.c -lm -o stelle.out #lo salvo come stelle.out
./stelle.out #lo eseguo

pkill -x gnuplot_qt,gnuplot

#!/usr/bin/gnuplot

gnuplot 'plot.plt'

rm stelle.out *.txt 
