#!/usr/bin/gnuplot

set term pngcairo enh color
set out "../images/specfun2.png"
file = '../examples/specfun2.txt'

load 'lines2.cfg'
load 'griddark.cfg'
set key top left

set xlabel "x"
set title "Spherical Harmonic ALFs"
plot file us 1:2 w li ti "Y_{21}", \
     file us 1:3 w li ti "dY_{21}/d{/Symbol \161}", \
     file us 1:4 w li ti "d^2 Y_{21}/d{/Symbol \161}^2"
