#!/bin/sh
#
# Run all tests under valgrind/memcheck
#
# First, compile the library with:
# > ./configure --disable-shared CFLAGS="-g -Wall"
# > make ; make check
# to compile the library and tests without optimization

outfile="valgrind.out"
rm -f $outfile

for testprog in $(ls */test); do
  echo "Running valgrind on ${testprog}"
  valgrind ${testprog} >> $outfile 2>&1
done

echo "output file is $outfile"
