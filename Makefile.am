## Process this file with automake to produce Makefile.in

SUBDIRS = utils err random randist fft roots specfunc statistics siman doc

include_HEADERS = gsl_complex.h gsl_math.h

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4"
