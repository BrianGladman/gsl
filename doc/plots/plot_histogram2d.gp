#!/usr/bin/env gnuplot

set term pngcairo enh col

unset key
set pointsize 0.3

set out "../images/histogram2d.png"
file = "../examples/histogram2d.txt"

set xrange [0:1]
set yrange [0:1]

set xlabel "x"
set ylabel "y"

set title "Distribution of simulated events"
plot file us 1:2 w p
