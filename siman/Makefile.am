## Process this file with automake to produce Makefile.in

CFLAGS=-g -Wall

check_PROGRAMS = siman_test siman_tsp
lib_LTLIBRARIES = libgslsiman.la

TESTS = siman_test_driver.sh
EXTRA_DIST = siman_test_driver.sh

siman_test_SOURCES = siman_test.c
siman_test_LDADD = libgslsiman.la ../random/libgslrandom.la ../utils/libutils.a
siman_test_DEPENDENCIES = libgslsiman.la

siman_tsp_SOURCES = siman_tsp.c
siman_tsp_LDADD = libgslsiman.la ../random/libgslrandom.la ../utils/libutils.a
siman_tsp_DEPENDENCIES = libgslsiman.la

CLEANFILES = siman_test.lout1 siman_test.lout2

libgslsiman_la_SOURCES = siman.c

include_HEADERS = gsl_siman.h

INCLUDES = -I$(srcdir) -I$(srcdir)/../random
check_INCLUDES = -I. -I../random
