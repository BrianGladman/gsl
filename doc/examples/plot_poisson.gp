#!/usr/bin/gnuplot

set term post eps enh color solid
set out "sparse_poisson.eps"

set xlabel "x"
set ylabel "u(x)"

load 'lines2.cfg'
load 'griddark.cfg'

set title "1D Poisson solution - offset deliberately added to GSL solution"
plot 'poisson.txt' us 1:($2 + 0.01) w li lw 4 ti "GSL + {/Symbol \145}", \
     'poisson.txt' us 1:3 w li lt 3 lw 4 ti "Analytic"
