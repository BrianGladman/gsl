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
 (cd gsl ; make LN_S=mv ; perl -i ../../msvc/mkdllheaders.pl *.h ; )
 (cd doc ; 
  test -e gsl-ref.info && rm -f *.info *.info-* ; 
  mkdir html ; 
  cd html ; 
  ../../../doc/texi2html -htmlhelp -verbose ../gsl-ref.texi ; )
  cp -a ../msvc . ;
 (cd msvc ; rm -rf usr *~ demo/*~ ; make ; cp -a gsl-ref.hhp ../doc/html; )
) 

if [ $? != 0 ] ; then
   echo failed to build;
   exit 1;
fi;


test -e $ver.zip && mv -b $ver.zip $ver.zip.old
zip -l -r $ver.zip $ver -x '*~' -x '*/CVS/*' -x '*/.*'
rm -r $ver


