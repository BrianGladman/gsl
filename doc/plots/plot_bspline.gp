#!/usr/bin/env gnuplot

set term pngcairo enh col

unset key
set pointsize 1

set out "../images/bspline.png"
file = '../examples/bspline.txt'

set xrange [0:15]
set xlabel "x"
set ylabel "y"

plot file index 0 us 1:2 w p lc rgb "black", \
     file index 1 us 1:2 w li lc rgb "red"
