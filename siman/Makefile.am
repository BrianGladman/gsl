## Process this file with automake to produce Makefile.in

CFLAGS=-g -Wall

check_PROGRAMS = siman_test siman_tsp
lib_LIBRARIES = libgsl_siman.a

siman_test_SOURCES = siman_test.c
siman_test_LDADD = libgsl_siman.a -L$(srcdir)/../random -lgslrandom
siman_test_DEPENDENCIES = libgsl_siman.a

siman_tsp_SOURCES = siman_tsp.c
siman_tsp_LDADD = libgsl_siman.a -L$(srcdir)/../random -lgslrandom
siman_tsp_DEPENDENCIES = libgsl_siman.a

libgsl_siman_a_SOURCES = siman.c

include_HEADERS = gsl_siman.h

INCLUDES = -I$(srcdir) -I$(srcdir)/../random
check_INCLUDES = -I. -I../random



