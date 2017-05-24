#!/usr/bin/env gnuplot

set term pngcairo enh mono

unset key
set pointsize 1

set out "../images/histogram.png"
file = "histogram.dat"

set xrange [-100:100]
set yzeroaxis

plot file us 1:3 w boxes
