#!/bin/sh
#
# Run all tests under valgrind/memcheck
#
# First, compile the library with:
# > make CFLAGS="-g -Wall" check
# to compile the library and tests without optimization

for testprog in `ls */test`; do
  echo "Running valgrind on ${testprog}"
  valgrind ${testprog}
done
