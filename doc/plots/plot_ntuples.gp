#!/usr/bin/env gnuplot

set term pngcairo enh mono

unset key

set out "../images/ntuple.png"
file = "../examples/ntuple.txt"

plot file us 1:3 w steps
