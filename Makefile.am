## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libgsl.a

DEP_FILES=""

SUBDIRS = err ieee-utils vector matrix blas specfunc rng randist fft roots statistics siman integration interpolation histogram monte doc

include_HEADERS = gsl_config.h gsl_complex.h gsl_math.h templates_on.h templates_off.h

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -Werror -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -Dinline='' -g -O4"

# Extend 'make all' through automake to build a single libgsl.a file

# The three following targets work around automake

libgsl.a.o:

libgsl.a.c:

.PHONY: libgsl.a.o libgsl.a.c

# The next target builds a new libgsl.a if any of the subdir/lib.a's
# are newer It has to be kept up to date manually, (C-u M-! echo
# */lib*.a) does the trick in emacs. There might be a problem if the
# build directory is different from the source directory.

libgsl.a: blas/libgslblas.a blas/libgslblascblas.a blas/libgslblasnative.a err/libgslerr.a err/libgsltest.a fft/libgslfft.a histogram/libgslhistogram.a ieee-utils/libgslieeeutils.a integration/libgslintegration.a interpolation/libgslinterpolation.a matrix/libgslmatrix.a monte/libgslmonte.a randist/libgslrandist.a rng/libgslrng.a roots/libgslroots.a siman/libgslsiman.a specfunc/libgslspecfunc.a statistics/libgslstatistics.a vector/libgslvector.a
	@rm -f libgsl.a
	for dir in $(SUBDIRS) ; \
	do (cd $$dir ; \
	  for lib in lib*.a ; \
	  do test -f $$lib \
	  && echo ar q libgsl.a `$(AR) t $$lib` \
	  && $(AR) q ../libgsl.a `$(AR) t $$lib` || true ; \
	  done ) ; \
	done 
	$(RANLIB) libgsl.a 
