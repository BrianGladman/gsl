lib_LTLIBRARIES = libgslerr.la libgsltest.la

include_HEADERS = gsl_errno.h gsl_message.h gsl_test.h

libgslerr_la_SOURCES = error.c stream.c message.c strerror.c warn.c
libgsltest_la_SOURCES = test.c

check_PROGRAMS = test_errnos

TESTS = test_errnos

test_errnos_LDADD = libgslerr.la libgsltest.la
