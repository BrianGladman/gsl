## Process this file with automake to produce Makefile.in

CFLAGS=-g -Wall

noinst_PROGRAMS = rtest
lib_LIBRARIES = gsl_roots

rtest_SOURCES = rtest.c
rtest_LDADD = libgsl_roots.a
rtest_DEPENDENCIES = libgsl_roots.a

gsl_roots_SOURCES = newton.c

include_HEADERS = roots.h
