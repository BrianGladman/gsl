noinst_LIBRARIES = libgslintegration.a

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err

libgslintegration_a_SOURCES = integration.c \
	qk15.c qk15.h \
	qk21.c qk21.h \
	qk31.c qk31.h \
	qk41.c qk41.h \
	qk51.c qk51.h \
	qk61.c qk61.h \
	qk.c \
	qng.c qng.h \
	qage.c qage_impl.c qpsrt.c qpsrt.h \
	qagse.c	qagse_impl.c qelg.c \
	err.c err.h max.h \
	workspace.c

include_HEADERS = gsl_integration.h

TESTS = test
check_PROGRAMS = test

test_SOURCES = test.c tests.c tests.h
test_LDADD = libgslintegration.a
test_LDADD = libgslintegration.a ../err/libgslerr.a ../err/libgsltest.a 
