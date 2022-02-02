# Example 5

set term pngcairo enh col size 1200,800
set out "../images/bspline_gram1.png"
file = '../examples/bspline_gram1.txt'
load 'lines.cfg'

# transparency factor
alpha = 0.3

set title "Fully regularized spline fit"
plot file index 1 us 1:($3-$4):($3+$4) w filledcu fs solid 0.3 lc rgb "red" ti "", \
     file index 1 us 1:($5-$6):($5+$6) w filledcu fs solid 0.3 lc rgb "blue" ti "", \
     file index 0 us 1:2 w p ps 1 pt 7 lc rgb "black" ti "Data", \
     file index 1 us 1:2 w li lw mylw lc rgb green_025 ti "Exact", \
     file index 1 us 1:3 w li lw mylw lc rgb "red" ti "Unregularized", \
     file index 1 us 1:5 w li lw mylw lc rgb "blue" ti "Regularized"

set term pngcairo enh col size 1200,800
set out "../images/bspline_gram2.png"
file = '../examples/bspline_gram2.txt'
load 'lines.cfg'

set title "Partially regularized spline fit"
plot file index 1 us 1:($3-$4):($3+$4) w filledcu fs solid 0.3 lc rgb "red" ti "", \
     file index 1 us 1:($5-$6):($5+$6) w filledcu fs solid 0.3 lc rgb "blue" ti "", \
     file index 0 us 1:2 w p ps 1 pt 7 lc rgb "black" ti "Data", \
     file index 1 us 1:2 w li lw mylw lc rgb green_025 ti "Exact", \
     file index 1 us 1:3 w li lw mylw lc rgb "red" ti "Unregularized", \
     file index 1 us 1:5 w li lw mylw lc rgb "blue" ti "Regularized"
