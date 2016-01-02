#!/usr/bin/gnuplot

set term post eps enh color solid

set out "nlfit2.eps"

set xlabel "x_1"
set ylabel "x_2"
set cblabel "{/Symbol \143}^2"

set pm3d map interp 20,20
set xrange [-1.2:1.2]

load 'lines2.cfg'
load 'moreland18.pal'

splot 'nlfit2.txt' index 0 us 1:2:3, \
      'nlfit2.txt' index 1 us 1:2:(0) w lp ps 1.5 pt 7 lt 1 ti "Without geodesic acceleration", \
      'nlfit2.txt' index 2 us 1:2:(0) w lp ps 1.5 pt 9 lt 2 ti "With geodesic acceleration"
