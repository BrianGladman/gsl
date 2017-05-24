#!/usr/bin/env gnuplot

set term pngcairo enh mono

unset key
set pointsize 1

set out "../images/dwt.png"
file = '../examples/dwt.txt'

set xzeroaxis
set xrange [0:256]
set xtics 0,32,256

set multiplot layout 2,1

plot file us 0:1 w li

plot file us 0:2 w li

unset multiplot
