#! /bin/bash

for CC in "gcc" "egcc"
do
for OPT in "-g -Wall" "-g -O2 -Wall"
do
make CC="$CC" CFLAGS="$OPT" test_sf
for mode in "" "double-precision;mask-all"
do
GSL_IEEE_MODE="$mode" ./test_sf
done
make clean
done
done
