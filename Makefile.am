## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libgsl.a

SUBDIRS = err ieee-utils vector matrix blas specfunc rng randist fft roots statistics siman integration interpolation histogram monte doc

include_HEADERS = gsl_config.h gsl_complex.h gsl_math.h templates_on.h templates_off.h

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -Werror -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -Dinline='' -g -O4"

# Extend 'make all' through automake to build a single libgsl.a file

# The three following targets work around automake

libgsl.a.o:

libgsl.a.c:

.PHONY: libgsl.a.o libgsl.a.c libgsl.a

# The next target builds a new libgsl.a if any of the subdir/lib.a's are newer

libgsl.a: 
	@for dir in $(SUBDIRS) ; do (cd $$dir ; for lib in lib*.a ; \
	  do test -f $$lib \
	  && test $$lib -nt ../libgsl.a \
	  && rm -f ../libgsl.a || true ; \
	  done ) ; \
	done	
	@test ! -f libgsl.a && \
	for dir in $(SUBDIRS) ; \
	do (cd $$dir ; \
	  for lib in lib*.a ; \
	  do test -f $$lib \
	  && echo ar q libgsl.a `ar t $$lib` \
	  && ar q ../libgsl.a `ar t $$lib` || true ; \
	  done ) ; \
	done || echo libgsl.a is up to date
