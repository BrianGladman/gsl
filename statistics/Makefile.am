## Process this file with automake to produce Makefile.in

#AUTOMAKE_OPTIONS=no-dependencies

CFLAGS=-g -Wall

lib_LIBRARIES = gsl_statistics
noinst_PROGRAMS = jtest

jtest_LDADD = libgsl_statistics.a
jtest_DEPENDENCIES = libgsl_statistics.a

gsl_statistics_SOURCES = mean.c ttest.c

jtest_SOURCES = jtest.c

include_HEADERS = mean.h
