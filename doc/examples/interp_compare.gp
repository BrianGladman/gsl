#!/usr/bin/gnuplot

set term post eps enh color solid

set out "interp_compare.eps"

set style line 1 pt 9 ps 1.5 lt -1
set style line 2 lt 1 lw 4
set style line 3 lt 2 lw 4
set style line 4 lt 3 lw 4

set key inside bottom right

plot 'interp_compare.txt' index 0 us 1:2 ls 1 ti "Data", \
     'interp_compare.txt' index 1 us 1:2 ls 2 wi li ti "Cubic", \
     'interp_compare.txt' index 1 us 1:3 ls 3 wi li ti "Akima", \
     'interp_compare.txt' index 1 us 1:4 ls 4 wi li ti "Steffen"
