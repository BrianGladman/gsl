#! /bin/bash
mkdir =gsl.build
cd =gsl.build
../configure --disable-shared
make clean
echo VERSION: `sh ./gsl-config --version` > typescript.make
echo running make
make -k CFLAGS="-g -O2 -Wall" >>typescript.make 2>&1 
echo running make check
make -k check CFLAGS="-g -O2 -Wall" >>typescript.make 2>&1 
../scripts/knownproblems.pl < typescript.make > ../KNOWN-PROBLEMS
