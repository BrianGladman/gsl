#!/bin/sh

./plot_bspline.gp
./plot_cheb.gp
./plot_dwt.gp
./plot_fft.gp
./plot_histogram.gp
./plot_histogram2d.gp
./plot_interp_compare.gp
./plot_lls.gp
./plot_min.gp
./plot_multimin.gp
./plot_nls.gp
./plot_ntuples.gp
./plot_ode.gp
./plot_qrng.gp
./plot_randist.gp
./plot_siman.gp

graph -T png < ../examples/interp.txt > ../images/interp.png
graph -T png < ../examples/interpp.txt > ../images/interpp.png

# LLS fitting example
for n in data fit hi lo ;
   do
     grep "^$n" ../examples/fitting.txt | cut -d: -f2 > $n ;
   done
graph -T png -X x -Y y -y 0 20 -m 0 -S 2 -Ie data -S 0 -I a -m 1 fit -m 2 hi -m 2 lo > ../images/fit-wlinear.png
