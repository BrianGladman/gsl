noinst_LIBRARIES = libgslintegration.a

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../ieee-utils -I$(srcdir)/../err

libgslintegration_a_SOURCES = \
	qk15.c qk15.h \
	qk21.c qk21.h \
	qk31.c qk31.h \
	qk41.c qk41.h \
	qk51.c qk51.h \
	qk61.c qk61.h \
	qk.c \
	qng.c qng.h \
	qag.c qag_impl.c qpsrt.c qpsrt.h \
	qags.c qags_impl.c qelg.c \
	qagi.c qagiu.c qagil.c \
	qagp.c qagp_impl.c \
	err.c err.h \
	workspace.c partition.c \
	qcheb.c qc25c.c

include_HEADERS = gsl_integration.h

TESTS = test
check_PROGRAMS = test

test_SOURCES = test.c tests.c tests.h
test_LDADD = libgslintegration.a
test_LDADD = libgslintegration.a ../ieee-utils/libgslieeeutils.a ../err/libgslerr.a ../err/libgsltest.a 
