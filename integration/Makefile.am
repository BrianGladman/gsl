pkglib_LIBRARIES = libgslintegration.a

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
	workspace.c \
	qcheb.c qawc.c \
	qmomo.c qaws.c \
	qmomof.c qawo.c \
	qawf.c 

pkginclude_HEADERS = gsl_integration.h
noinst_HEADERS = integration.h qpsrt.c qpsrt2.c qelg.c qc25c.c qc25s.c qc25f.c ptsort.c util.c 


TESTS = test
check_PROGRAMS = test

test_SOURCES = test.c tests.c tests.h
test_LDADD = libgslintegration.a ../ieee-utils/libgslieeeutils.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
