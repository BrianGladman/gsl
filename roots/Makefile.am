## Process this file with automake to produce Makefile.in

CFLAGS=-g -Wall

check_PROGRAMS = rtest
lib_LIBRARIES = libgsl_roots.a

rtest_SOURCES = rtest.c
rtest_LDADD = libgsl_roots.a
rtest_DEPENDENCIES = libgsl_roots.a

libgsl_roots_a_SOURCES = newton.c

include_HEADERS = gsl_roots.h
