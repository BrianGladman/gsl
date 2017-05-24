#!/usr/bin/env gnuplot

set term pngcairo enh mono size 640,320

set out "../images/min-interval.png"

unset key
set pointsize 1

set xrange [-3:3]
set yrange [0:*]
f(x) = x**2 + 3
str1 = sprintf("<echo -1.5 %f", f(-1.5))
str2 = sprintf("<echo 1.1 %f", f(1.1))
str3 = sprintf("<echo 2.2 %f", f(2.2))

set label 1 "(a)" at -1.5,4.5 center
set label 2 "(x)" at 1.1,3.2 center
set label 3 "(b)" at 2.2,7.0 center
plot f(x) w li, \
     str1 w p pt 7, \
     str2 w p pt 7, \
     str3 w p pt 7
unset label 1
unset label 2
unset label 3
