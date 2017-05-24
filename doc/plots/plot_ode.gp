#!/usr/bin/env gnuplot

set term pngcairo enh mono

unset key

set out "../images/ode-vdp.png"
file = "../examples/ode-initval.txt"

set yrange [-5:5]

plot file us 1:2 w li
