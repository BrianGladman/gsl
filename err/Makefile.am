noinst_LIBRARIES = libgslerr.a libgsltest.a

include_HEADERS = gsl_errno.h gsl_message.h gsl_test.h

libgslerr_a_SOURCES = errno.c
libgsltest_a_SOURCES = test.c

check_PROGRAMS = test_errnos

TESTS = test_errnos

test_errnos_LDADD = libgslerr.a libgsltest.a


