#!/usr/bin/env gnuplot

set term pngcairo enh mono size 600,800

unset key
set yrange [-0.5:2]
set xzeroaxis
set pointsize 1

set out "../images/fft-complex-radix2.png"
file = "../examples/fft.txt"

set multiplot layout 2,1

plot file index 0 us 1:2 w p pt 7
plot file index 1 us 1:2 w p pt 7

unset multiplot

# Second plot for real mixed-radix FFT

set term pngcairo enh mono size 640,480
set yrange [*:*]
set pointsize 0.5
set out "../images/fft-real-mixedradix.png"
file = "../examples/fftreal.txt"

plot file us 1:2 w lp pt 7
