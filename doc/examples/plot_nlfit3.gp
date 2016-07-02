#!/usr/bin/gnuplot

set term post eps enh color solid

set out "nlfit3.eps"

set xlabel "x_1"
set ylabel "x_2"

unset surface
set contour
set cntrparam levels 50
set palette maxcolors 0
load 'jet.pal'

set table 'cntrs.dat'
splot 'nlfit3.txt' index 0 us 1:2:3 w li
unset table

load 'lines2.cfg'
set view map
set key bottom right tc variable font "Arial-Bold" opaque
set xrange [-5:15]
set yrange [-5:15]

set title "Minimizing the Branin function"
plot 'cntrs.dat' us 1:2:3 w li palette ti "", \
     'nlfit3.txt' index 1 us 1:2 w lp lw 6 ps 1.5 pt 7 lt 1 ti "LM", \
     'nlfit3.txt' index 2 us 1:2 w lp lw 6 ps 1.5 pt 9 lt 2 ti "LM + geodesic acceleration", \
     'nlfit3.txt' index 3 us 1:2 w lp lw 6 ps 1.5 pt 6 lt 3 ti "Dogleg", \
     'nlfit3.txt' index 4 us 1:2 w lp lw 6 ps 1.5 pt 2 lt 4 ti "Double Dogleg", \
     'nlfit3.txt' index 5 us 1:2 w lp lw 6 ps 1.5 pt 3 lt 5 ti "2D Subspace"
