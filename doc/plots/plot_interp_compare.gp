#!/usr/bin/gnuplot

set term pngcairo enh color solid

set out "../images/interp_compare.png"
file = "../examples/interp_compare.txt"

set style line 1 pt 9 ps 1.5 lt -1
set style line 2 lt 1 lw 4
set style line 3 lt 2 lw 4
set style line 4 lt 3 lw 4

load 'grid.cfg'

set key inside bottom right

plot file index 1 us 1:2 ls 2 wi li ti "Cubic", \
     file index 1 us 1:3 ls 3 wi li ti "Akima", \
     file index 1 us 1:4 ls 4 wi li ti "Steffen", \
     file index 0 us 1:2 ls 1 ti "Data"
