noinst_LIBRARIES = libgslintegration.a

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../ieee-utils -I$(srcdir)/../err

libgslintegration_a_SOURCES = \
	qk15.c qk21.c qk31.c qk41.c qk51.c qk61.c \
	qk.c \
	qng.c qng.h \
	qag.c qag_impl.c \
	qags.c qags_impl.c \
	qagi.c qagiu.c qagil.c \
	qagp.c qagp_impl.c \
	err.c err.h \
	workspace.c

#	qcheb.c qc25c.c qawc.c qawc_impl.c

include_HEADERS = gsl_integration.h
noinst_HEADERS = integration.h qpsrt.c qelg.c

TESTS = test
check_PROGRAMS = test

test_SOURCES = test.c tests.c tests.h
test_LDADD = libgslintegration.a
test_LDADD = libgslintegration.a ../ieee-utils/libgslieeeutils.a ../err/libgslerr.a ../err/libgsltest.a 
