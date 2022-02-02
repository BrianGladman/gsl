# Example 6

reset

set term pngcairo enh col size 1200,800
set out "../images/bspline_lsend.png"
file = '../examples/bspline_lsend.txt'
load 'lines.cfg'
set grid

# transparency factor
alpha = 0.3

inset_y = 0.47

set multiplot

set format y "%.1f"

set object 1 ellipse center -0.98,0.05 size 0.1,0.1
set object 2 ellipse center 0.98,0.05 size 0.1,0.1
set arrow 1 from -0.97,0.12 to screen 0.15,inset_y front
set arrow 2 from 0.97,0.12 to screen 0.92,inset_y front

set grid
set title "Runge function spline fit"
plot file index 1 us 1:($3-$4):($3+$4) w filledcu fs solid 0.3 lc rgb "red" ti "", \
     file index 1 us 1:($5-$6):($5+$6) w filledcu fs solid 0.3 lc rgb "blue" ti "", \
     file index 0 us 1:2 w p ps 1 pt 7 lc rgb "black" ti "Data", \
     file index 1 us 1:2 w li lw mylw lc rgb green_025 ti "Exact", \
     file index 1 us 1:3 w li lw mylw lc rgb "red" ti "Unregularized", \
     file index 1 us 1:5 w li lw mylw lc rgb "blue" ti "Regularized"

unset object
unset arrow
unset grid
unset key
unset title
set xtics 0.05
set format y "%.2f"

set origin 0.07,inset_y
set size 0.35, 0.35
plot [-1:-0.9] file index 1 us 1:($3-$4):($3+$4) w filledcu fs solid 0.3 lc rgb "red" ti "", \
     file index 1 us 1:($5-$6):($5+$6) w filledcu fs solid 0.3 lc rgb "blue" ti "", \
     file index 0 us 1:2 w p ps 1 pt 7 lc rgb "black" ti "Data", \
     file index 1 us 1:2 w li lw mylw lc rgb green_025 ti "Exact", \
     file index 1 us 1:3 w li lw mylw lc rgb "red" ti "Unregularized", \
     file index 1 us 1:5 w li lw mylw lc rgb "blue" ti "Regularized"

set origin 0.62,inset_y
set size 0.35, 0.35
plot [0.9:1] file index 1 us 1:($3-$4):($3+$4) w filledcu fs solid 0.3 lc rgb "red" ti "", \
     file index 1 us 1:($5-$6):($5+$6) w filledcu fs solid 0.3 lc rgb "blue" ti "", \
     file index 0 us 1:2 w p ps 1 pt 7 lc rgb "black" ti "Data", \
     file index 1 us 1:2 w li lw mylw lc rgb green_025 ti "Exact", \
     file index 1 us 1:3 w li lw mylw lc rgb "red" ti "Unregularized", \
     file index 1 us 1:5 w li lw mylw lc rgb "blue" ti "Regularized"

unset multiplot
