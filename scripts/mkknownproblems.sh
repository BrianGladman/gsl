#! /bin/bash
DIR=${DIR:-/tmp/gsl.build}
SRCDIR=`pwd`/..			# assume we build in $(srcdir)/scripts
echo building in $DIR 
mkdir $DIR
cd $DIR
#../configure --disable-shared
#../configure
${SRCDIR}/configure
make clean
echo VERSION: `sh ./gsl-config --version` > typescript.make
echo Running make >>typescript.make 2>&1 
make -k CFLAGS="-g -O2 -Wall" >>typescript.make 2>&1 
echo Running make check >>typescript.make 2>&1 
make -k check CFLAGS="-g -O2 -Wall" >>typescript.make 2>&1 
echo Running make check in double-precision >>typescript.make 2>&1 
export GSL_IEEE_MODE=double-precision,mask-all
make -k check CFLAGS="-g -O2 -Wall" >>typescript.make 2>&1 
${SRCDIR}/scripts/knownproblems.pl < typescript.make > ${SRCDIR}/KNOWN-PROBLEMS
