#!/usr/bin/gnuplot

set term post eps enh color solid size 7,3.5

set out "multilarge.eps"

set xlabel "t"
set ylabel "f(t)"

set multiplot layout 1,2

set yrange [-1:*]

set title "Not regularized ({/Symbol \154} = 0)"
plot 'largefit.out' us 1:2 w li lw 4 ti "Exact", \
     'largefit.out' us 1:3 w li lw 4 ti "TSQR", \
     'largefit.out' us 1:4 w li lw 4 ti "Normal"

unset key

set title "Regularized ({/Symbol \154} = 1)"
plot 'largefit2.out' us 1:2 w li lw 4 ti "Exact", \
     'largefit2.out' us 1:3 w li lw 4 ti "TSQR", \
     'largefit2.out' us 1:4 w li lw 4 ti "Normal"

unset multiplot
