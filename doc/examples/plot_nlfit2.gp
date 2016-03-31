#!/usr/bin/gnuplot

set term post eps enh color solid

set out "nlfit2.eps"

set xlabel "x_1"
set ylabel "x_2"

unset surface
set contour
set cntrparam levels 50
set palette maxcolors 0
load 'jet.pal'

set table 'cntrs.dat'
splot 'nlfit2.txt' index 0 us 1:2:3 w li
unset table

load 'lines2.cfg'
set view map
set key tc variable font "Arial-Bold" opaque
set xrange [-1.2:1.2]

plot 'cntrs.dat' us 1:2:3 w li palette ti "", \
     'nlfit2.txt' index 1 us 1:2 w lp lw 4 ps 1.5 pt 7 lt 1 ti "LM", \
     'nlfit2.txt' index 2 us 1:2 w lp lw 4 ps 1.5 pt 9 lt 3 ti "LM + geodesic acceleration"
