#!/usr/bin/env gnuplot

set term pngcairo enh col size 1600,800

load 'grid.cfg'

set xlabel "time"
set point 1

set out "../images/movstat1.png"
file = "../examples/movstat1.txt"
load 'lines2.cfg'

mylw = 4
set key bottom right

plot file us 1:2 w lp lc rgb "gray" lw mylw ti "Original signal", \
     file us 1:3 w li lw mylw ti "Moving mean", \
     file us 1:4 w li lw mylw ti "Moving minimum", \
     file us 1:5 w li lw mylw ti "Moving maximum"

set out "../images/movstat2.png"
file = "../examples/movstat2.txt"

load 'lines.cfg'
set key top left

set multiplot layout 2,1

plot file us 1:2 w li lc rgb "black" ti "x(t)"

plot file us 1:3 w li ti "True sigma", \
     file us 1:4 w li ti "MAD", \
     file us 1:5 w li ti "IQR", \
     file us 1:6 w li ti "S_n", \
     file us 1:7 w li ti "Q_n", \
     file us 1:8 w li lc rgb "gray" ti "sigma"

unset multiplot

set out "../images/movstat3.png"
file = "../examples/movstat3.txt"

plot file us 0:1 w li lw 3 lc rgb "black" ti "Data", \
     file us 0:2 w li lw mylw lt 7 ti "Trimmed mean"

unset multiplot
