noinst_LTLIBRARIES = libgslfit.la 

pkginclude_HEADERS = gsl_fit.h

INCLUDES= -I$(top_builddir)

libgslfit_la_SOURCES = linear.c

check_PROGRAMS = test

#noinst_PROGRAMS = demo
#demo_LDADD = libgslfit.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la

TESTS = test

test_SOURCES = test.c
test_LDADD = libgslfit.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la
test_DEPENDENCIES = libgslfit.la
