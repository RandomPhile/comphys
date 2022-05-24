#!/bin/bash

gcc errore.c -lm -o errore.out #lo salvo come stelle.out
./errore.out #lo eseguo

pkill -x gnuplot_qt,gnuplot

#!/usr/bin/gnuplot

gnuplot 'plot_errore.plt'

rm errore.out *.txt 
