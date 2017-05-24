#!/usr/bin/env gnuplot

set term pngcairo enh mono

set key inside right bottom
set pointsize 1

set out "../images/cheb.png"
file = '../examples/cheb.txt'

plot file us 1:2 w li lw 2 ti "Original function", \
     file us 1:3 w li lw 2 ti "10th order", \
     file us 1:4 w li lw 2 ti "40th order"
