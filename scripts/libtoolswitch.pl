#!/usr/bin/perl -i

@files = ("configure.in","Makefile.am",<*/Makefile.am>) ;

$switch = $ARGV[0] ;

@ARGV = @files ;

if ($switch eq 'on') {
    print "switching libtool on...\n" ;
    while (<>) {
        s/pkglib_pkglib/pkglib/g ;
	#s/^pkglib_LIBRARIES/pkglib_LTLIBRARIES/g ;
	s/^pkglib_LIBRARIES/noinst_LTLIBRARIES/g ;
	s/^\#?\s*libgsl_a_LIBADD/libgsl_la_LIBADD/g ;
	s/^\#?\s*lib_LTLIBRARIES = libgsl.la/lib_LTLIBRARIES = libgsl.la/g ;
	s/^\#?\s*libgsl_a_SOURCES/libgsl_la_SOURCES/g ;
	s/lib(\w+)\.a/lib$1.la/g ;
	s/lib(\w+)_a/lib$1_la/g ;
#	s/libutils\.la/libutils.a/g ;  # keep libutils as .a always
#	s/libutils_la/libutils_a/g ;
	s/(\w+)\.o/$1.lo/g ;
	s/^AC_PROG_RANLIB/#AC_PROG_RANLIB/ ;
	s/^\#AM_PROG_LIBTOOL/AM_PROG_LIBTOOL/ ;
	print ;
    }
} elsif ($switch eq 'off') {
    print "switching libtool off...\n" ;
    while (<>) {
        s/^(pkglib_LTLIBRARIES = libgsl.la)/\# \1/;
        s/^(lib_LTLIBRARIES = libgsl.la)/\# \1/;
	s/^pkglib_LTLIBRARIES/pkglib_LIBRARIES/g ;
	s/^noinst_LTLIBRARIES/pkglib_LIBRARIES/g ;
	s/^lib_LTLIBRARIES = libgsl.la/\# lib_LTLIBRARIES = libgsl.la/g ;
	s/^libgsl_la_SOURCES/\# libgsl_a_SOURCES/g ;
	s/libgsl_la_LIBADD/\# libgsl_a_LIBADD/g ;
	s/lib(\w+)\.la/lib$1.a/g ;
	s/lib(\w+)_la/lib$1_a/g ;
#	s/libutils\.la/libutils.a/g ; # keep libutils as .a always
#	s/libutils_la/libutils_a/g ;
	s/(\w+)\.lo/$1.o/g ;
	s/^\#AC_PROG_RANLIB/AC_PROG_RANLIB/ ;
	s/^AM_PROG_LIBTOOL/\#AM_PROG_LIBTOOL/ ;
	print ;
    }
} else {
    print <<EOF ;
libtoolswitch.pl [on|off]

Use this command in the top-level directory of gsl
It modifies the following files to use libtool or non-libtool style macros:

EOF
    print join(", ",@files),"\n" ;
}

