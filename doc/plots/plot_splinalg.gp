#!/usr/bin/gnuplot

set term pngcairo enh color
set out "../images/sparse_poisson.png"
file = '../examples/poisson.txt'

set xlabel "x"
set ylabel "u(x)"

load 'lines2.cfg'
load 'griddark.cfg'

set title "1D Poisson solution - offset deliberately added to GSL solution"
plot file us 1:($2 + 0.02) w li lw 4 ti "GSL + {/Symbol \145}", \
     file us 1:3 w li lt 3 lw 4 ti "Analytic"
