## Process this file with automake to produce Makefile.in

#AUTOMAKE_OPTIONS=no-dependencies

#CFLAGS=-g -Wall

noinst_PROGRAMS = jtest
lib_LIBRARIES = libgslstatistics.a

jtest_SOURCES = jtest.c
jtest_LDADD = libgslstatistics.a
jtest_DEPENDENCIES = libgslstatistics.a

libgslstatistics_a_SOURCES = mean.c variance.c test.c

#include_HEADERS = mean.h test.h variance.h
include_HEADERS = gsl_bstats.h
