## Process this file with automake to produce Makefile.in

#AUTOMAKE_OPTIONS=no-dependencies

CFLAGS=-g -Wall

noinst_PROGRAMS = jtest
lib_LIBRARIES = libgsl_statistics.a

jtest_SOURCES = jtest.c
jtest_LDADD = libgsl_statistics.a
jtest_DEPENDENCIES = libgsl_statistics.a

libgsl_statistics_a_SOURCES = mean.c variance.c test.c

#include_HEADERS = mean.h test.h variance.h
include_HEADERS = gsl_bstats.h
