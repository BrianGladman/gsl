## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libgsl.a

SUBDIRS = err ieee-utils vector matrix blas specfunc rng randist fft roots statistics siman integration interpolation histogram monte doc

include_HEADERS = gsl_config.h gsl_complex.h gsl_math.h templates_on.h templates_off.h

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -Werror -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -Dinline='' -g -O4"

# Extend 'make all' through automake to build a single libgsl.a file

libgsl.a.o:

libgsl.a.c:

.PHONY: libgsl.a.o libgsl.a.c libgsl.a

libgsl.a: 
	-test -f libgsl.a && touch -r libgsl.a stamp-h1
	for dir in $(SUBDIRS) ; \
	do (cd $$dir ; \
	  for lib in lib*.a ; \
	  do test -f $$lib \
	  && (test ! -f ../stamp-h1 || test $$lib -nt ../stamp-h1)\
	  && echo ar q libgsl.a `ar t $$lib` \
	  && ar q ../libgsl.a `ar t $$lib` || true ; \
	  done ) ; \
	done	
	for dir in $(SUBDIRS) ; do (cd $$dir ; for lib in lib*.a ; \
	  do test -f $$lib && test $$lib -nt ../libgsl.a \
	  && touch -r ../libgsl.a $$lib || true ; \
	  done ) ; \
	done	
	-rm -f stamp-h1
