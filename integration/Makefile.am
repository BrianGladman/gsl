lib_LIBRARIES = libgslintegration.a

libgslintegration_a_SOURCES = integration.c \
	qk15.c qk21.c qk31.c qk41.c qk51.c qk61.c qk.c

include_HEADERS = gsl_integration.h
check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslintegration.a
