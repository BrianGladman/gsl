#!/usr/bin/env gnuplot

set term pngcairo enh mono

unset key
set pointsize 1

set out "../images/multimin.png"
file = '../examples/multimin.txt'

p0 = 1
p1 = 2
p2 = 10
p3 = 20
p4 = 30
f(x,y) = p2 * (x - p0)**2 + p3 * (y - p1)**2 + p4

set xrange [0:8]
set yrange [0:8]

set contour
set view map
unset surface
set cntrparam levels 15
set isosamples 100

set table $datatable
splot f(x,y)
unset table

plot $datatable w p pt 7 ps 0.1, \
     file us 2:3 w lp pt 7
