lib_LIBRARIES = libgslintegration.a

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
	qage.c qpsrt.c \
	err.c err.h

include_HEADERS = gsl_integration.h
check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslintegration.a
test_LDADD = libgslintegration.a ../err/libgslerr.a ../err/libgsltest.a 
