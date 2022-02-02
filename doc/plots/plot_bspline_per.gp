# Example 7

reset

set term pngcairo enh col size 1200,800
set out "../images/bspline_per.png"
file = '../examples/bspline_per.txt'
load 'lines3.cfg'
set grid
mylw = 6

set title "Non-periodic and periodic spline fits"
plot file index 0 us 1:2 w p ps 1 pt 7 lc rgb "black" ti "Data", \
     file index 1 us 1:2 w li lt 6 lw mylw ti "Non-periodic", \
     file index 1 us 1:3 w li lt 4 lw mylw ti "Periodic"
