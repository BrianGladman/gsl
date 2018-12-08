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
unset xlabel

## Gaussian example 1

set term pngcairo enh col size 1200,1000

set out "../images/gaussfilt.png"
file = "../examples/gaussfilt.txt"

set key inside bottom right

set multiplot layout 2,1

plot file index 0 us 0:1 w li lw mylw ti "Gaussian kernel for alpha = 0.5", \
     file index 0 us 0:2 w li lw mylw ti "Gaussian kernel for alpha = 3", \
     file index 0 us 0:3 w li lw mylw ti "Gaussian kernel for alpha = 10"

plot file index 1 us 0:1 w li lw mylw lc rgb "gray" ti "Data", \
     file index 1 us 0:2 w li lw mylw lt 1 ti "Smoothed data for alpha = 0.5", \
     file index 1 us 0:3 w li lw mylw lt 2 ti "Smoothed data for alpha = 3", \
     file index 1 us 0:4 w li lw mylw lt 3 ti "Smoothed data for alpha = 10"

unset multiplot

## Gaussian example 2

set term pngcairo enh col size 1000,1200

set out "../images/gaussfilt2.png"
file = "../examples/gaussfilt2.txt"

set key inside top left

set multiplot layout 4,1

plot file us 0:1 w li lc rgb "black" ti "Signal", \
     file us 0:2 w li lt 3 lw mylw ti "Gaussian smoothed signal"

plot file us 0:3 w li lc rgb "black" ti "First differenced signal"

plot file us 0:4 w li lc rgb "black" ti "First order Gaussian smoothed signal"

plot file us 0:5 w li lc rgb "black" ti "Second order Gaussian smoothed signal"

unset multiplot
