CC=gcc
COPT=-g
CFLAGS=$(COPT)

sig2lamp: sig2lamp.o erfc.o
	$(CC) -o $@ $@.o erfc.o -lm

t-erf: t-erf.o erfc.o
	$(CC) -o $@ $@.o erfc.o -lm


