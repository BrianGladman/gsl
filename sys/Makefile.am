noinst_LTLIBRARIES = libgslsys.la 

libgslsys_la_SOURCES = minmax.c prec.c hypot.c log1p.c expm1.c

INCLUDES = -I$(top_builddir)

check_PROGRAMS = test
TESTS = test
test_sources = test.c
test_LDADD = libgslsys.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la
