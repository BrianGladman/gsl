## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libgslstatistics.a
include_HEADERS = gsl_statistics.h

noinst_PROGRAMS = test

TESTS = test

test_SOURCES = test.c test.h
test_LDADD = libgslstatistics.a ../err/libgsltest.a
test_DEPENDENCIES = libgslstatistics.a

libgslstatistics_a_SOURCES = \
mean.c variance.c \
mean.h variance.h \
minmax.c
