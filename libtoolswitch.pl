#!/usr/bin/perl -i

@files = ("configure.in","Makefile.am",<*/Makefile.am>) ;

$switch = $ARGV[0] ;

@ARGV = @files ;

if ($switch eq 'on') {
    print "switching libtool on...\n" ;
    while (<>) {
	s/lib_LIBRARIES/lib_LTLIBRARIES/g ;
	s/lib(\w+)\.a/lib$1.la/g ;
	s/lib(\w+)_a/lib$1_la/g ;
	s/libutils\.la/libutils.a/g ;  # keep libutils as .a always
	s/libutils_la/libutils_a/g ;
	s/(\w+)\.o/$1.lo/g ;
	s/^AC_PROG_RANLIB/#AC_PROG_RANLIB/ ;
	s/^\#AM_PROG_LIBTOOL/AM_PROG_LIBTOOL/ ;
	print ;
    }
} else {
    print "switching libtool off...\n" ;
    while (<>) {
	s/lib_LTLIBRARIES/lib_LIBRARIES/g ;
	s/lib(\w+)\.la/lib$1.a/g ;
	s/lib(\w+)_la/lib$1_a/g ;
	s/libutils\.la/libutils.a/g ; # keep libutils as .a always
	s/libutils_la/libutils_a/g ;
	s/(\w+)\.lo/$1.o/g ;
	s/^\#AC_PROG_RANLIB/AC_PROG_RANLIB/ ;
	s/^AM_PROG_LIBTOOL/\#AM_PROG_LIBTOOL/ ;
	print ;
    }
}    

