#!/usr/bin/gnuplot

set term post eps enh color solid size 8,4

set out "regularized2.eps"

set style line 1 pt 9 ps 1.5 lt -1
set style line 2 lt 1 lw 4
set style line 3 lt 2 lw 4
set style line 4 lt 3 lw 4

set logscale x
set logscale y
set format x "10^{%L}"
set format y "10^{%L}"
unset key

set multiplot layout 1,2

set xlabel "residual norm ||y - X c||"
set ylabel "solution norm ||c||"

set xrange [2.3:3.2]
set xtics ("2.3" 2.3, "3.2" 3.2)

set title "L-curve"
plot 'fitreg2.txt' index 0 us 2:3 w lp lt -1 pt 7, \
     'fitreg2.txt' index 1 us 1:2 w p ps 3 pt 6 lt 7 lw 2

set xlabel "{/Symbol \154}"
set ylabel "G({/Symbol \154})"
set xrange [*:*]
set xtics auto

set title "GCV curve"
plot 'fitreg2.txt' index 0 us 1:4 w lp lt -1 pt 7, \
     'fitreg2.txt' index 2 us 1:2 w p ps 3 pt 6 lt 7 lw 2

unset multiplot
