## Process this file with automake to produce Makefile.in
CFLAGS=-g -Wall
# CC=gcc

noinst_PROGRAMS = jtest
noinst_LIBRARIES = gsl_statistics

LDADD = libgsl_statistics.a -lm

gsl_statistics_SOURCES = mean.c

jtest_SOURCES = jtest.c

include_HEADERS = mean.h
