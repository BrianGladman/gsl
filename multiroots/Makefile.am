# -*-makefile-*-

pkglib_LIBRARIES = libgslmultiroots.a 

pkginclude_HEADERS = gsl_multiroots.h

noinst_HEADERS = multiroots.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err

libgslroots_a_SOURCES = 

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c 
test_LDADD = libgslmultiroots.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_DEPENDENCIES = libgslmultiroots.a
