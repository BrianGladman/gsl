## Process this file with automake to produce Makefile.in

SUBDIRS = err ieee-utils specfunc rng vector matrix randist fft roots statistics siman integration interpolation histogram monte blas doc

include_HEADERS = gsl_config.h gsl_complex.h gsl_math.h templates_on.h templates_off.h

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -Werror -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -Dinline='' -g -O4"

# a flaky attempt to build a single libgsl.a file (can someone improve on it?)

libgsl.a: 
	for dir in $(SUBDIRS) ;\
	do \
	(cd $$dir ; \
	for lib in lib*.a ; \
	do \
	test -f $$lib && ar r ../libgsl.a `ar t $$lib` || true ; \
	done ) ;\
	done

# flaky attempt to generalize the above flaky version to work with libtool.
# Not yet, baby. 

#libgsl.a: 
#	for dir in $(SUBDIRS) ;\
#	do \
#	(cd $$dir ; \
#	for lib in ./.libs/lib*.a ; \
#	do \
#	test -f $$lib && ar r ../libgsl.a `ar t $$lib|egrep -v "^__\.SYMDEF" ` || true ; \
#	done ) ;\
#	done
