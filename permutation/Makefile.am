pkglib_LIBRARIES = libgslpermutation.a 

check_PROGRAMS = test

pkginclude_HEADERS = gsl_permutation.h

INCLUDES= -I$(top_builddir)

TESTS = test

test_LDADD = ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a

test_DEPENDENCIES = libgslpermutation.a

test_SOURCES = test.c

CLEANFILES = test.txt test.dat

libgslpermutation_a_SOURCES = init.c file.c permutation.c
