#!/usr/bin/env gnuplot

set term pngcairo enh mono

unset key
set pointsize 1

set out "../images/qrng.png"
file = "../examples/qrng.txt"

plot file us 1:2 w p pt 7
