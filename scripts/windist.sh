#! /bin/bash

set -e  # quit on error
set -x  # display commands

tarfile=$1

if [ ! -f "$tarfile" ] ; then
   echo $1 must be a tar file ;
fi;

ver=${tarfile%%.tar.gz}

tar xvfz $tarfile
( cd $ver;
  ./configure;
 (cd gsl ; make)
 (cd doc ; mkdir html ; 
  cd html ; 
  ../../../doc/texi2html -htmlhelp -verbose ../gsl-ref.texi ; )
 cp -a ../msvc .
 (cd msvc ; make ; cp -a gsl-ref.hhp ../doc/html; )
)

zip -l -r $ver.zip $ver
rm -r $ver


