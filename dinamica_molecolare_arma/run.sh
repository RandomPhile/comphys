export CPATH=/opt/homebrew/include
export LIBRARY_PATH=/opt/homebrew/lib

#!/bin/bash
pkill -x gnuplot, gnuplot_qt
g++ main.cpp -o main -std=c++11 -O2 -larmadillo -o main.out
./main.out
time ./out/a.out
