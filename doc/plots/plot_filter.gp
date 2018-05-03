#!/usr/bin/env gnuplot

set term pngcairo enh col size 1600,800

load 'grid.cfg'

set xlabel "time"
set point 1
mylw = 4

set out "../images/impulse.png"
file = "../examples/impulse.txt"

plot file us 1:2 w li ti "Data", \
     file us 1:3 w li ti "Filtered data", \
     file us 1:4 w li ti "Upper limit", \
     file us 1:5 w li ti "Lower limit", \
     file us 1:($6 == 1 ? $2 : 1/0) w p ps 1.5 pt 4 ti "Outliers"

set out "../images/filt_edge.png"
file = "../examples/filt_edge.txt"

set xlabel "time (s)"
plot file us 1:2 w li lw mylw lc rgb "gray" ti "Data", \
     file us 1:3 w li lw mylw ti "Standard Median Filter", \
     file us 1:4 w li lw mylw ti "Recursive Median Filter"
