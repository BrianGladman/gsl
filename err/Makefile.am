pkglib_LIBRARIES = libgslerr.a libgsltest.a

pkginclude_HEADERS = gsl_errno.h gsl_message.h gsl_test.h

libgslerr_a_SOURCES = error.c stream.c message.c strerror.c warn.c
libgsltest_a_SOURCES = test_results.c

check_PROGRAMS = test_errnos

TESTS = test_errnos

test_errnos_LDADD = libgslerr.a libgsltest.a ../utils/libutils.a
