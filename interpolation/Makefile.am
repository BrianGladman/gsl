lib_LIBRARIES = libgslinterpolation.a 

# check_PROGRAMS = test_interp

include_HEADERS = 		\
	gsl_interp.h

noinst_HEADERS = 		\
	bsearch.h		\
	tridiag.h


INCLUDES = -I$(full_top_srcdir) -I$(full_top_srcdir)/err

#TESTS = test_interp

#test_sf_LDADD = libgslinterpolation.a ../err/libgslerr.a ../err/libgsltest.a
#test_sf_DEPENDENCIES = libgslinterpolation.a
#test_sf_SOURCES = test_interp.c

libgslinterpolation_a_SOURCES = 	\
	accel.c				\
	bsearch.c			\
	cspline.c			\
	interp.c			\
	linear.c			\
	tridiag.c


strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
