## Process this file with automake to produce Makefile.in

#AUTOMAKE_OPTIONS=no-dependencies

CFLAGS=-g -Wall

noinst_PROGRAMS = jtest
lib_LIBRARIES = gsl_statistics

jtest_SOURCES = jtest.c
jtest_LDADD = libgsl_statistics.a
jtest_DEPENDENCIES = libgsl_statistics.a

gsl_statistics_SOURCES = mean.c variance.c test.c

include_HEADERS = mean.h test.h variance.h
