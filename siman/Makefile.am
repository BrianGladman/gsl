## Process this file with automake to produce Makefile.in

CFLAGS=-g -Wall

check_PROGRAMS = siman_test siman_tsp
lib_LIBRARIES = libgslsiman.a

TESTS = siman_test_driver.sh
EXTRA_DIST = siman_test_driver.sh

siman_test_SOURCES = siman_test.c
siman_test_LDADD = libgslsiman.a -L../random -lgslrandom -L../utils -lutils
siman_test_DEPENDENCIES = libgslsiman.a

siman_tsp_SOURCES = siman_tsp.c
siman_tsp_LDADD = libgslsiman.a -L../random -lgslrandom -L../utils  -lutils
siman_tsp_DEPENDENCIES = libgslsiman.a

CLEANFILES = siman_test.out1 siman_test.out2

libgslsiman_a_SOURCES = siman.c

include_HEADERS = gsl_siman.h

INCLUDES = -I$(srcdir) -I$(srcdir)/../random
check_INCLUDES = -I. -I../random
