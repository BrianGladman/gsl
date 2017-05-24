#!/usr/bin/gnuplot
#
# Plot the output of interp2d.c with gnuplot

set term post eps color enh
set output "interp2d.eps"

set pm3d map

# similar to matlab color palette
set palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee',\
                      4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000',\
                      8 '#7f0000')

set xlabel "x"
set ylabel "y"
set label 1 "z = 0" at graph 0.01,0.05 left front tc rgb "white"
set label 2 "z = 1" at graph 0.01,0.95 left front tc rgb "white"
set label 3 "z = 1" at graph 0.99,0.05 right front tc rgb "white"
set label 4 "z = 0.5" at graph 0.99,0.95 right front tc rgb "white"
splot 'interp2d.txt' us 1:2:3
