## Process this file with automake to produce Makefile.in

SUBDIRS = utils err ieee-utils random vector matrix randist fft roots specfunc statistics siman integration histogram doc

include_HEADERS = gsl_complex.h gsl_math.h

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -Werror -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -Dinline='' -g -O4"

# a flaky attempt to build a single libgsl.la file (can someone improve on it?)

libgsl.la: 
	for dir in $(SUBDIRS) ;\
	do \
	(cd $$dir ; \
	for lib in ./.libs/lib*.a ; \
	do \
	test -f $$lib && ar r ../libgsl.la `ar t $$lib|egrep -v "^__\.SYMDEF" ` || true ; \
	done ) ;\
	done

