pkglib_LIBRARIES = libgslinterpolation.a 

check_PROGRAMS = test

pkginclude_HEADERS = 		\
	gsl_interp.h

noinst_HEADERS = 		\
	bsearch.h		\
	tridiag.h


INCLUDES = -I$(full_top_srcdir) -I$(full_top_srcdir)/err

TESTS = test

test_LDADD = libgslinterpolation.a ../err/libgslerr.a ../err/libgsltest.a
test_DEPENDENCIES = libgslinterpolation.a
test_SOURCES = test.c

libgslinterpolation_a_SOURCES = 	\
	accel.c				\
	akima.c				\
	bsearch.c			\
	cspline.c			\
	interp.c			\
	linear.c			\
	tridiag.c
