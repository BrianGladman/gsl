noinst_LTLIBRARIES = libgsldht.la 

pkginclude_HEADERS = gsl_dht.h

INCLUDES= -I$(top_builddir)

TESTS = test

check_PROGRAMS = test

test_LDADD = libgsldht.la ../specfunc/libgslspecfunc.la ../sys/libgslsys.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la

test_DEPENDENCIES = libgsldht.la
test_SOURCES = test_dht.c

libgsldht_la_SOURCES = dht_transform.c
