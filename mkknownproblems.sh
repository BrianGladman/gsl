#! /bin/bash
make clean
make -k CFLAGS="-g -O2 -Wall" 2>&1 | tee typescript.make
make -k check CFLAGS="-g -O2 -Wall" 2>&1 | tee -a typescript.make
./knownproblems.pl < typescript.make > KNOWN-PROBLEMS
