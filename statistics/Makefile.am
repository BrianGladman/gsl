## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libgslstatistics.a

include_HEADERS = gsl_statistics.h

libgslstatistics_a_SOURCES = \
mean.c variance.c \
mean.h variance.h

TESTS = test

test_SOURCES = test.c test.h
test_LDADD = libgslstatistics.a
test_DEPENDENCIES = libgslstatistics.a

