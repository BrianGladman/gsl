## Process this file with automake to produce Makefile.in

CFLAGS=-g -Wall

check_PROGRAMS = siman_test
lib_LIBRARIES = libgsl_siman.a

siman_test_SOURCES = siman_test.c
siman_test_LDADD = libgsl_siman.a -lgslrandom
siman_test_DEPENDENCIES = libgsl_siman.a

libgsl_siman_a_SOURCES = siman.c

include_HEADERS = gsl_siman.h

INCLUDES = -I. -I../random
