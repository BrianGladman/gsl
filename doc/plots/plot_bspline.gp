#!/usr/bin/env gnuplot

# Example 1

set term pngcairo enh col size 2000,1200

file_spline = '../examples/bspline_knots.txt'
file_knots = '../examples/bspline_knots_data.txt'

load 'grid.cfg'
unset key
set xtics 0.1

set ylabel "B(x)"
set format x ""

set out "../images/bspline_knots.png"

set multiplot layout 2,2 columnsfirst margins 0.05,0.98,0.05,0.95 spacing 0.05,0.05

# Plot linear B-splines

idx = 0
set title "Linear B-splines"
plot for [i=2:12] file_spline index idx us 1:(column(i)) w li, \
     file_knots index idx us 1:(0.0) w p pt 6 ps 3 lc rgb "black"

# Plot quadratic B-splines

set format x "%g"
set xlabel "x"
idx = 1
set title "Quadratic B-splines"
plot for [i=2:13] file_spline index idx us 1:(column(i)) w li, \
     file_knots index idx us 1:(0.0) w p pt 6 ps 3 lc rgb "black"

unset xlabel
set format x ""

# Plot cubic B-splines

set label 1 at 0.02, 0.8 'B_1(x)' front tc lt 1
set label 2 at 0.03, 0.63 'B_2(x)' front tc lt 2
set label 3 at 0.09, 0.63 'B_3(x)' front tc lt 3
set label 4 at 0.18, 0.7 'B_4(x)' front tc lt 4
set label 5 at 0.28, 0.7 'B_5(x)' front tc lt 5
set label 6 at 0.38, 0.7 'B_6(x)' front tc lt 6
set label 7 at 0.48, 0.7 'B_7(x)' front tc lt 7
set label 8 at 0.58, 0.7 'B_8(x)' front tc lt 8
set label 9 at 0.68, 0.7 'B_9(x)' front tc lt 9
set label 10 at 0.78, 0.7 'B_{10}(x)' front tc lt 10
set label 11 at 0.86, 0.63 'B_{11}(x)' front tc lt 11
set label 12 at 0.93, 0.63 'B_{12}(x)' front tc lt 12
set label 13 at 0.93, 0.8 'B_{13}(x)' front tc lt 13

idx = 2
set title "Cubic B-splines"
plot for [i=2:14] file_spline index idx us 1:(column(i)) w li, \
     file_knots index idx us 1:(0.0) w p pt 6 ps 3 lc rgb "black"

unset label

set format x "%g"
set xlabel "x"
idx = 3
set title "Quartic B-splines"
plot for [i=2:15] file_spline index idx us 1:(column(i)) w li, \
     file_knots index idx us 1:(0.0) w p pt 6 ps 3 lc rgb "black"

unset label

unset multiplot
set xtics auto
unset ylabel
unset title

# Example 2

set term pngcairo enh col size 2000,1200

file = '../examples/bspline_deriv.txt'

load 'grid.cfg'
unset key
set xtics 0.1

set out "../images/bspline_deriv.png"

set multiplot layout 2,2 columnsfirst margins 0.05,0.98,0.05,0.95 spacing 0.05,0.05

# Plot cubic B-splines

set label 1 at 0.02, 0.8 'B_1(x)' front tc lt 1
set label 2 at 0.07, 0.63 'B_2(x)' front tc lt 2
set label 3 at 0.20, 0.63 'B_3(x)' front tc lt 3
set label 4 at 0.38, 0.7 'B_4(x)' front tc lt 4
set label 5 at 0.58, 0.7 'B_5(x)' front tc lt 5
set label 6 at 0.76, 0.63 'B_6(x)' front tc lt 6
set label 7 at 0.89, 0.63 'B_7(x)' front tc lt 7
set label 8 at 0.93, 0.8 'B_8(x)' front tc lt 8

idx = 1
set ylabel "B(x)"
set format x ""
unset xlabel
set title "Cubic B-splines"
plot for [i=2:14] file index idx us 1:(column(i)) w li, \
     file index 0 us 1:(0.0) w p pt 6 ps 3 lc rgb "black"

unset label

idx = 2
set format x "%g"
set xlabel "x"
set ylabel "d/dx B(x)"
set title "Cubic B-splines first derivatives"
plot for [i=2:14] file index idx us 1:(column(i)) w li, \
     file index 0 us 1:(0.0) w p pt 6 ps 3 lc rgb "black"

idx = 3
set format x ""
unset xlabel
set ylabel "d^2/dx^2 B(x)"
set title "Cubic B-splines second derivatives"
plot for [i=2:14] file index idx us 1:(column(i)) w li, \
     file index 0 us 1:(0.0) w p pt 6 ps 3 lc rgb "black"

idx = 4
set format x "%g"
set xlabel "x"
set ylabel "d^3/dx^3 B(x)"
set title "Cubic B-splines third derivatives"
plot for [i=2:14] file index idx us 1:(column(i)) w li, \
     file index 0 us 1:(0.0) w p pt 6 ps 3 lc rgb "black"

unset multiplot
set xtics auto
unset ylabel
unset title

# Example 3

set term pngcairo enh col size 1200,800
set out "../images/bspline_lsbreak.png"

set key top right inside
set xrange [0:15]
set xtics 1
file = '../examples/bspline_lsbreak.txt'
set title "Cubic spline least squares fit"
plot file index 0 us 1:2 w p ps 1 lc rgb "black" ti "Data", \
     file index 1 us 1:2 w li lw 6 lt 3 ti "40 breakpoints", \
     file index 1 us 1:3 w li lw 6 lt 5 ti "10 breakpoints"

# Example 4

set term pngcairo enh col size 1200,800
set out "../images/bspline_lsorder.png"

set key top right inside
set xrange [0:15]
set xtics 1
mylw = 4
file = '../examples/bspline_lsorder.txt'
set title "Various order spline fits"
plot file index 0 us 1:2 w p ps 1 lc rgb "black" ti "Data", \
     file index 1 us 1:2 w li lw mylw ti "Order 1", \
     file index 1 us 1:3 w li lw mylw ti "Order 2", \
     file index 1 us 1:4 w li lw mylw ti "Order 3", \
     file index 1 us 1:5 w li lw mylw ti "Order 4", \
     file index 1 us 1:6 w li lw mylw ti "Order 5"

set xrange [*:*]
set xtics auto

load 'plot_bspline_gram.gp'
load 'plot_bspline_lsend.gp'
load 'plot_bspline_per.gp'

# Example 8

set term pngcairo enh col size 1200,800
set out "../images/bspline_proj.png"
file = '../examples/bspline_proj.txt'
load 'lines3.cfg'
mylw = 6

set title "Projection onto B-spline basis"
plot file us 1:2 w li lt 6 lw mylw ti "Exact", \
     file us 1:($3+0.1) w li lt 4 lw mylw ti "Projection"

load 'plot_bspline_interp.gp'
