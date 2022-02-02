# Example 9

reset

set term pngcairo enh col size 1200,800
set out "../images/bspline_interp.png"
file = '../examples/bspline_interp.txt'
load 'lines3.cfg'
set grid
mylw = 6

set title "Interpolated spline"
plot file index 1 us 1:2 w li lt 1 lw mylw ti "Spline", \
     file index 0 us 1:2 w p ps 2 pt 7 lc rgb "black" ti "Data"
