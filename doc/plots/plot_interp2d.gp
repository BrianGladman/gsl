#!/usr/bin/gnuplot
#
# Plot the output of interp2d.c with gnuplot

set term pngcairo color enh
set output "../images/interp2d.png"

file = '../examples/interp2d.txt'

set pm3d map

load 'jet.pal'

set xlabel "x"
set ylabel "y"
set label 1 "z = 0" at graph 0.01,0.05 left front tc rgb "white"
set label 2 "z = 1" at graph 0.01,0.95 left front tc rgb "white"
set label 3 "z = 1" at graph 0.99,0.05 right front tc rgb "white"
set label 4 "z = 0.5" at graph 0.99,0.95 right front tc rgb "white"
splot file us 1:2:3
