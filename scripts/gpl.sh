#! /bin/bash
PATH=./scripts:$PATH

addgpl.pl "Gerard Jungman" gsl_machine.h gsl_mode.h gsl_precision.h
addgpl.pl "Brian Gough" version.c gsl_version.h

addgpl.pl "Gerard Jungman, Brian Gough" gsl_math.h complex/gsl_complex.h
addgpl.pl "Jorma Olavi Tähtinen, Brian Gough" complex/gsl_complex_math.h complex/math.c

addgpl.pl "Jim Davies, Brian Gough" statistics/*.{c,h}
addgpl.pl "Mark Galassi" siman/*.{c,h}
addgpl.pl "Michael Booth" monte/*.{c,h}
addgpl.pl "Fabrice Rossi" multimin/*.{c,h}
addgpl.pl "Reid Priedhorsky, Brian Gough" roots/*.{c,h}
addgpl.pl "Thomas Walter, Brian Gough" sort/*.{c,h}

for dir in blas eigen specfunc dht interpolation ode-initval
do
  addgpl.pl "Gerard Jungman" $dir/*.{c,h}
done

for dir in sys err block vector matrix linalg sum
do
  addgpl.pl "Gerard Jungman, Brian Gough" $dir/*.{c,h}
done

for dir in rng randist
do
  addgpl.pl "James Theiler, Brian Gough" $dir/*.{c,h}
done

addgpl.pl "Tim Mooney" ieee-utils/fp-{tru64,irix,aix}.c

for dir in fft integration min multiroots poly histogram permutation ieee-utils
do
  addgpl.pl "Brian Gough" $dir/*.{c,h}
done


