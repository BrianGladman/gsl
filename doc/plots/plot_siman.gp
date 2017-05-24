#!/usr/bin/env gnuplot

set term pngcairo enh mono

unset key
set pointsize 1

file = "../examples/siman.txt"

set out "../images/siman-test.png"
set yrange [1.34:1.4]
set xlabel "generation"
set ylabel "position"
plot file us 1:4 w li

set out "../images/siman-energy.png"
set yrange [-0.88:-0.83]
set xlabel "generation"
set ylabel "energy"
plot file us 1:5 w li

file = "../examples/siman_tsp.txt"

set out "../images/siman-12-cities.png"
set yrange [3300:6500]
set xlabel "generation"
set ylabel "distance"
set title "TSP - 12 southwest cities"
plot file us 1:18 w li

set out "../images/siman-initial-route.png"
set title "TSP - initial order"
set yrange [*:*]
set xlabel "longitude (- means west)"
set ylabel "latitude"
plot "<grep initial_city_coord ".file." | awk '{print $2,$3}'" us 1:2 w lp pt 7

set out "../images/siman-final-route.png"
set title "TSP - final order"
set yrange [*:*]
set xlabel "longitude (- means west)"
set ylabel "latitude"
plot "<grep final_city_coord ".file." | awk '{print $2,$3}'" us 1:2 w lp pt 7
