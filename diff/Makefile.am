noinst_LTLIBRARIES = libgsldiff.la

INCLUDES= -I$(top_builddir)

libgsldiff_la_SOURCES = diff.c
	
pkginclude_HEADERS = gsl_diff.h

TESTS = test
check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = \
	   libgsldiff.la \
	   ../multimin/libgslmultimin.la \
	   ../vector/libgslvector.la \
	   ../block/libgslblock.la \
	   ../err/libgslerr.la \
	   ../test/libgsltest.la
