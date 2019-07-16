#!/usr/bin/env gnuplot

set term pngcairo enh mono

unset key
set pointsize 1

## fitting 2 example

# best fit curve
f(x) = 1.02318 + 0.956201*x + 0.876796*x**2

set out "../images/fit-wlinear2.png"
plot [0:2] f(x), \
           '../examples/fitting2.txt' us 1:2:4 w errorbars pt 7

## regularized fit example 1

set term pngcairo enh col size 1000,480

set out "../images/regularized.png"
file = '../examples/fitreg.txt'

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

set title "L-curve"
plot file index 0 us 2:3 w lp lt -1 pt 7, \
     file index 1 us 1:2 w p ps 3 pt 6 lt 7 lw 2

set xlabel "{/Symbol \154}"
set ylabel "G({/Symbol \154})"
set yrange [1e-4:*]

set title "GCV curve"
plot file index 0 us 1:4 w lp lt -1 pt 7, \
     file index 2 us 1:2 w p ps 3 pt 6 lt 7 lw 2

unset multiplot

## regularized fit example 2

set out "../images/regularized2.png"
file = '../examples/fitreg2.txt'

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
plot file index 0 us 2:3 w lp lt -1 pt 7, \
     file index 1 us 1:2 w p ps 3 pt 6 lt 7 lw 2

set xlabel "{/Symbol \154}"
set ylabel "G({/Symbol \154})"
set xrange [*:*]
set xtics auto

set title "GCV curve"
plot file index 0 us 1:4 w lp lt -1 pt 7, \
     file index 2 us 1:2 w p ps 3 pt 6 lt 7 lw 2

unset multiplot

unset logscale x
unset logscale y
set yrange [*:*]
set xrange [*:*]
set format x "%g"
set format y "%g"
unset title
unset xlabel
unset ylabel

## robust fitting example

set term pngcairo enh col size 640,480

set out "../images/robust.png"
file = '../examples/robfit.txt'

set key left inside top

plot file us 1:2 w p pt 9 ti "Data", \
     file us 1:3 w li lw 4 lt 3 ti "Robust", \
     file us 1:4 w li lw 4 lt 4 ti "OLS"

## multilarge example

set term pngcairo enh col size 1300,1000

set out "../images/multilarge.png"
file = '../examples/largefit.txt'
file2 = '../examples/largefit2.txt'

set xlabel "t"
set ylabel "f(t)"

set multiplot layout 2,2 rowsfirst

load 'lines2.cfg'
set yrange [0:4]
set key top left inside

set title "Not regularized ({/Symbol \154} = 0)"
plot file us 1:2 index 0 w p pt 7 ps 0.5 lc rgb "black" ti "Data", \
     file us 1:2 index 3 w li lw 6 ti "Exact", \
     file us 1:3 index 3 w li lw 6 ti "TSQR", \
     file us 1:4 index 3 w li lw 4 ti "Normal"

unset key

set title "Regularized ({/Symbol \154} = 1e-5)"
plot file2 us 1:2 index 0 w p pt 7 ps 0.5 lc rgb "black" ti "Data", \
     file2 us 1:2 index 3 w li lw 6 ti "Exact", \
     file2 us 1:3 index 3 w li lw 6 ti "TSQR", \
     file2 us 1:4 index 3 w li lw 4 ti "Normal"

set xrange [1e1:1e3]
set yrange [1e-1:1e10]
load 'xylogon.cfg'
load 'lines.cfg'
set xlabel "residual norm ||y - X c||"
set ylabel "solution norm ||c||"

set title "L-curve computed from TSQR method"
plot file us 2:3 index 1 w li lw 6

set title "L-curve computed from normal equations method"
plot file us 2:3 index 2 w li lw 6

unset multiplot
