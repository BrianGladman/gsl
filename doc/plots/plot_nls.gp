#!/usr/bin/env gnuplot

set term pngcairo enh color

set pointsize 1

## exponential fitting example

# best fit curve
A = 4.79653
lambda = 1.43937
b = 1.00368
f(x) = A * exp(-lambda*x) + b

load 'grid.cfg'

set xlabel "time"
set out "../images/fit-exp.png"
plot '../examples/nlfit.txt' us 2:3:4 ls 1 w errorbars ti "Data", \
     f(x) ls 3 lw 4 ti "Model"

## geodesic acceleration example 1

set term pngcairo enh color size 800,600
set out "../images/nlfit2.png"
file = "../examples/nlfit2.txt"

set xlabel "x_1"
set ylabel "x_2"

unset surface
set contour
set cntrparam levels 50
set palette maxcolors 0
load 'jet.pal'

set table 'cntrs.dat'
splot file index 0 us 1:2:3 w li
unset table

load 'lines2.cfg'
set view map
set key tc variable font "Arial-Bold" opaque
set xrange [-1.2:1.2]

plot 'cntrs.dat' us 1:2:3 w li palette ti "", \
     file index 1 us 1:2 w lp lw 4 ps 1.5 pt 7 lt 1 ti "LM", \
     file index 2 us 1:2 w lp lw 4 ps 1.5 pt 9 lt 3 ti "LM + geodesic acceleration"

set xrange [*:*]

## geodesic acceleration example 2

set out "../images/nlfit2b.png"
file = "../examples/nlfit2b.txt"

load 'grid.cfg'

set xlabel "t"
set ylabel "y"
set point 1

plot file us 1:2 w p pt 7 lc "#000000" ti "Data", \
     file us 1:3 w li lt 5 lw 6 ti "Model"

## comparing TRS example

set out "../images/nlfit3.png"
file = "../examples/nlfit3.txt"

set xlabel "x_1"
set ylabel "x_2"

unset surface
set contour
set cntrparam levels 50
set palette maxcolors 0
load 'jet.pal'

set table 'cntrs.dat'
splot file index 0 us 1:2:3 w li
unset table

load 'lines2.cfg'
set view map
set key bottom right tc variable font "Arial-Bold" opaque
set xrange [-5:15]
set yrange [-5:15]

set title "Minimizing the Branin function"
plot 'cntrs.dat' us 1:2:3 w li palette ti "", \
     file index 1 us 1:2 w lp lw 6 ps 1.5 pt 7 lt 1 ti "LM", \
     file index 2 us 1:2 w lp lw 6 ps 1.5 pt 9 lt 2 ti "LM + geodesic acceleration", \
     file index 3 us 1:2 w lp lw 6 ps 1.5 pt 6 lt 3 ti "Dogleg", \
     file index 4 us 1:2 w lp lw 6 ps 1.5 pt 2 lt 4 ti "Double Dogleg", \
     file index 5 us 1:2 w lp lw 6 ps 1.5 pt 3 lt 5 ti "2D Subspace"
