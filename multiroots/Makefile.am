# -*-makefile-*-

pkglib_LIBRARIES = libgslmultiroots.a 

pkginclude_HEADERS = gsl_multiroots.h

#noinst_HEADERS = multiroots.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../matrix -I$(srcdir)/../vector -I$(srcdir)/../block  -I$(srcdir)/../err

libgslmultiroots_a_SOURCES = fdjac.c fsolver.c fdfsolver.c

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c 
test_LDADD = libgslmultiroots.a ../matrix/libgslmatrix.a ../vector/libgslvector.a ../block/libgslblock.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_DEPENDENCIES = libgslmultiroots.a
