#! /bin/bash
make clean
echo VERSION: `sh ./gsl-config --version` > typescript.make
make -k CFLAGS="-g -O2 -Wall" 2>&1 | tee -a typescript.make
make -k check CFLAGS="-g -O2 -Wall" 2>&1 | tee -a typescript.make
./scripts/knownproblems.pl < typescript.make > KNOWN-PROBLEMS
