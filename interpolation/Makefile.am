pkglib_LIBRARIES = libgslinterpolation.a 

check_PROGRAMS = test

pkginclude_HEADERS = 		\
	gsl_interp.h

noinst_HEADERS = 		\
	bsearch.h


INCLUDES = \
	-I$(full_top_srcdir)		\
	-I$(full_top_srcdir)/err 	\
	-I$(full_top_srcdir)/linalg	\
	-I$(srcdir)/../permutation      \
	-I$(full_top_srcdir)/vector	\
	-I$(full_top_srcdir)/block	\
	-I$(full_top_srcdir)/matrix

TESTS = test

test_LDADD = libgslinterpolation.a 	\
	../linalg/libgsllinalg.a 	\
	../vector/libgslvector.a	\
	../block/libgslblock.a		\
	../err/libgslerr.a 		\
	../err/libgsltest.a 		\
	../utils/libutils.a

test_DEPENDENCIES = libgslinterpolation.a
test_SOURCES = test.c

libgslinterpolation_a_SOURCES = 	\
	accel.c				\
	akima.c				\
	bsearch.c			\
	cspline.c			\
	interp.c			\
	linear.c			\
        integ_eval_macro.h
