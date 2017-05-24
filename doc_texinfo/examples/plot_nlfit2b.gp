#!/usr/bin/gnuplot

set term post eps enh color solid
set out "nlfit2b.eps"

load 'grid.cfg'

set xlabel "t"
set ylabel "y"
set point 1

plot 'nlfit2b.txt' us 1:2 w p pt 7 lc "#000000" ti "Data", \
     'nlfit2b.txt' us 1:3 w li lt 5 lw 6 ti "Model"
