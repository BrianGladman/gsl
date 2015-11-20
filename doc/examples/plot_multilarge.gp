#!/usr/bin/gnuplot

set term post eps enh color solid size 7,6

set out "multilarge.eps"

set xlabel "t"
set ylabel "f(t)"

set multiplot layout 2,2 rowsfirst

load 'lines2.cfg'
set yrange [0:4]
set key top left inside

set title "Not regularized ({/Symbol \154} = 0)"
plot 'largefit.txt' us 1:2 index 0 w p pt 7 ps 0.5 lc rgb "black" ti "Data", \
     'largefit.txt' us 1:2 index 3 w li lw 6 ti "Exact", \
     'largefit.txt' us 1:3 index 3 w li lw 6 ti "TSQR", \
     'largefit.txt' us 1:4 index 3 w li lw 6 ti "Normal"

unset key

set title "Regularized ({/Symbol \154} = 1e-6)"
plot 'largefit2.txt' us 1:2 index 0 w p pt 7 ps 0.5 lc rgb "black" ti "Data", \
     'largefit2.txt' us 1:2 index 3 w li lw 6 ti "Exact", \
     'largefit2.txt' us 1:3 index 3 w li lw 6 ti "TSQR", \
     'largefit2.txt' us 1:4 index 3 w li lw 6 ti "Normal"

set xrange [1e1:1e3]
set yrange [1e-1:1e10]
load 'xylogon.cfg'
load 'lines.cfg'
set xlabel "residual norm ||y - X c||"
set ylabel "solution norm ||c||"

set title "L-curve computed from TSQR method"
plot 'largefit.txt' us 2:3 index 1 w li lw 6

set title "L-curve computed from normal equations method"
plot 'largefit.txt' us 2:3 index 2 w li lw 6

unset multiplot
