## $Id$
## This file is a template for Makefile.am (which is a template for 
## Makefile.in, which is a template for Makefile...) 
## You can bootstrap your way to a regular
## Makefile with 'make -f Makefile.ami Makefile.am' followed by the usual
## automake/configure routines

## Include here the list of algorithms for random number generation
## The first one should be the "favorite" (ie, default)

ALLRAN= taus mrg cmrg uni uni32 zuf rand
FAVORITERAN= taus

include_HEADERS= gsl_ran.h gsl_ran_switch.h

gauss_SOURCES= rungauss.c gauss.c
test_gauss_SOURCES= test-gauss.c gauss.c
bench_gauss_SOURCES= bench-gauss.c gauss.c

gauss_LDADD=       $(FAVORITERAN).o $(FAVORITERAN)-gen.o
test_gauss_LDADD=  $(FAVORITERAN).o $(FAVORITERAN)-gen.o
bench_gauss_LDADD= $(FAVORITERAN).o $(FAVORITERAN)-gen.o

test_switch_SOURCES= test-switch.c switch.c 

switch.c gsl_ran_switch.h: switch.xc makeswitch
	makeswitch $(ALLRAN) 

## The user doesn't really _need_ these, but 'make all' does,
## and probably it is fair for the user to have them, since that
## permits him/her to contribute to the project!!
##
EXTRA_DIST= switch.xc makeswitch xxx-gen.xc xxx-state.xc \
	makeam Makefile.ami

CLEANFILES= tmp.out tmp.in

## The rest of the Makefile.am is generated from the script 'makeam'
## which includes numerous loops over $(ALLRAN)

Makefile.am: Makefile.ami makeam
	cat Makefile.ami    > Makefile.am
	./makeam $(ALLRAN) >> Makefile.am


## need to automake'ify...somehow

test-shuffle: test-shuffle.o $(FAVORITERAN).o $(FAVORITERAN)-gen.o
	$(LINK) -o $@ test-shuffle.o $(FAVORITERAN).o $(FAVORITERAN)-gen.o -lm

test-shuffle.o: shuffle.c
	$(COMPILE) -c -DMAIN -o $@ shuffle.c





bin_PROGRAMS= mrg test-mrg bench-mrg cmrg test-cmrg bench-cmrg taus test-taus bench-taus uni test-uni bench-uni uni32 test-uni32 bench-uni32 zuf test-zuf bench-zuf rand test-rand bench-rand gauss test-gauss bench-gauss test-switch 

lib_LIBRARIES= gslrandom gslranmrg gslrancmrg gslrantaus gslranuni gslranuni32 gslranzuf gslranrand 

gslranmrg_SOURCES= mrg.c mrg-gen.c gauss.c 
gslrancmrg_SOURCES= cmrg.c cmrg-gen.c gauss.c 
gslrantaus_SOURCES= taus.c taus-gen.c gauss.c 
gslranuni_SOURCES= uni.c uni-gen.c gauss.c 
gslranuni32_SOURCES= uni32.c uni32-gen.c gauss.c 
gslranzuf_SOURCES= zuf.c zuf-gen.c gauss.c 
gslranrand_SOURCES= rand.c rand-gen.c gauss.c 

gslrandom_SOURCES= switch.c shuffle.c mrg.c cmrg.c taus.c uni.c uni32.c zuf.c rand.c gauss.c 

test_switch_DEPENDENCIES= gsl_ran_switch.h mrg.o cmrg.o taus.o uni.o uni32.o zuf.o rand.o 
test_switch_LDADD= mrg.o cmrg.o taus.o uni.o uni32.o zuf.o rand.o 
mrg_SOURCES=             ran.c mrg.c
mrg_DEPENDENCIES=        mrg-state.c mrg.h mrg-gen.o
mrg_LDADD=               mrg-gen.o
test_mrg_SOURCES=        testran.c
test_mrg_DEPENDENCIES=   mrg.o mrg-gen.o mrg.h
test_mrg_LDADD=          mrg.o mrg-gen.o
bench_mrg_SOURCES=       benchran.c
bench_mrg_DEPENDENCIES=  mrg.o mrg-gen.o mrg.h
bench_mrg_LDADD=         mrg.o mrg-gen.o
cmrg_SOURCES=             ran.c cmrg.c
cmrg_DEPENDENCIES=        cmrg-state.c cmrg.h cmrg-gen.o
cmrg_LDADD=               cmrg-gen.o
test_cmrg_SOURCES=        testran.c
test_cmrg_DEPENDENCIES=   cmrg.o cmrg-gen.o cmrg.h
test_cmrg_LDADD=          cmrg.o cmrg-gen.o
bench_cmrg_SOURCES=       benchran.c
bench_cmrg_DEPENDENCIES=  cmrg.o cmrg-gen.o cmrg.h
bench_cmrg_LDADD=         cmrg.o cmrg-gen.o
taus_SOURCES=             ran.c taus.c
taus_DEPENDENCIES=        taus-state.c taus.h taus-gen.o
taus_LDADD=               taus-gen.o
test_taus_SOURCES=        testran.c
test_taus_DEPENDENCIES=   taus.o taus-gen.o taus.h
test_taus_LDADD=          taus.o taus-gen.o
bench_taus_SOURCES=       benchran.c
bench_taus_DEPENDENCIES=  taus.o taus-gen.o taus.h
bench_taus_LDADD=         taus.o taus-gen.o
uni_SOURCES=             ran.c uni.c
uni_DEPENDENCIES=        uni-state.c uni.h uni-gen.o
uni_LDADD=               uni-gen.o
test_uni_SOURCES=        testran.c
test_uni_DEPENDENCIES=   uni.o uni-gen.o uni.h
test_uni_LDADD=          uni.o uni-gen.o
bench_uni_SOURCES=       benchran.c
bench_uni_DEPENDENCIES=  uni.o uni-gen.o uni.h
bench_uni_LDADD=         uni.o uni-gen.o
uni32_SOURCES=             ran.c uni32.c
uni32_DEPENDENCIES=        uni32-state.c uni32.h uni32-gen.o
uni32_LDADD=               uni32-gen.o
test_uni32_SOURCES=        testran.c
test_uni32_DEPENDENCIES=   uni32.o uni32-gen.o uni32.h
test_uni32_LDADD=          uni32.o uni32-gen.o
bench_uni32_SOURCES=       benchran.c
bench_uni32_DEPENDENCIES=  uni32.o uni32-gen.o uni32.h
bench_uni32_LDADD=         uni32.o uni32-gen.o
zuf_SOURCES=             ran.c zuf.c
zuf_DEPENDENCIES=        zuf-state.c zuf.h zuf-gen.o
zuf_LDADD=               zuf-gen.o
test_zuf_SOURCES=        testran.c
test_zuf_DEPENDENCIES=   zuf.o zuf-gen.o zuf.h
test_zuf_LDADD=          zuf.o zuf-gen.o
bench_zuf_SOURCES=       benchran.c
bench_zuf_DEPENDENCIES=  zuf.o zuf-gen.o zuf.h
bench_zuf_LDADD=         zuf.o zuf-gen.o
rand_SOURCES=             ran.c rand.c
rand_DEPENDENCIES=        rand-state.c rand.h rand-gen.o
rand_LDADD=               rand-gen.o
test_rand_SOURCES=        testran.c
test_rand_DEPENDENCIES=   rand.o rand-gen.o rand.h
test_rand_LDADD=          rand.o rand-gen.o
bench_rand_SOURCES=       benchran.c
bench_rand_DEPENDENCIES=  rand.o rand-gen.o rand.h
bench_rand_LDADD=         rand.o rand-gen.o

mrg-state.c: xxx-state.xc
	sed s/gsl_ran_xxx_/gsl_ran_mrg_/g xxx-state.xc > $@
mrg.h: gsl_ran.h
	sed s/gsl_ran_/gsl_ran_mrg_/g gsl_ran.h > $@
mrg-gen.c: xxx-gen.xc
	sed s/xxx/mrg/g xxx-gen.xc > $@
cmrg-state.c: xxx-state.xc
	sed s/gsl_ran_xxx_/gsl_ran_cmrg_/g xxx-state.xc > $@
cmrg.h: gsl_ran.h
	sed s/gsl_ran_/gsl_ran_cmrg_/g gsl_ran.h > $@
cmrg-gen.c: xxx-gen.xc
	sed s/xxx/cmrg/g xxx-gen.xc > $@
taus-state.c: xxx-state.xc
	sed s/gsl_ran_xxx_/gsl_ran_taus_/g xxx-state.xc > $@
taus.h: gsl_ran.h
	sed s/gsl_ran_/gsl_ran_taus_/g gsl_ran.h > $@
taus-gen.c: xxx-gen.xc
	sed s/xxx/taus/g xxx-gen.xc > $@
uni-state.c: xxx-state.xc
	sed s/gsl_ran_xxx_/gsl_ran_uni_/g xxx-state.xc > $@
uni.h: gsl_ran.h
	sed s/gsl_ran_/gsl_ran_uni_/g gsl_ran.h > $@
uni-gen.c: xxx-gen.xc
	sed s/xxx/uni/g xxx-gen.xc > $@
uni32-state.c: xxx-state.xc
	sed s/gsl_ran_xxx_/gsl_ran_uni32_/g xxx-state.xc > $@
uni32.h: gsl_ran.h
	sed s/gsl_ran_/gsl_ran_uni32_/g gsl_ran.h > $@
uni32-gen.c: xxx-gen.xc
	sed s/xxx/uni32/g xxx-gen.xc > $@
zuf-state.c: xxx-state.xc
	sed s/gsl_ran_xxx_/gsl_ran_zuf_/g xxx-state.xc > $@
zuf.h: gsl_ran.h
	sed s/gsl_ran_/gsl_ran_zuf_/g gsl_ran.h > $@
zuf-gen.c: xxx-gen.xc
	sed s/xxx/zuf/g xxx-gen.xc > $@
rand-state.c: xxx-state.xc
	sed s/gsl_ran_xxx_/gsl_ran_rand_/g xxx-state.xc > $@
rand.h: gsl_ran.h
	sed s/gsl_ran_/gsl_ran_rand_/g gsl_ran.h > $@
rand-gen.c: xxx-gen.xc
	sed s/xxx/rand/g xxx-gen.xc > $@

BUILT_SOURCES= switch.c gsl_ran_switch.h mrg-state.c mrg.h mrg-gen.c cmrg-state.c cmrg.h cmrg-gen.c taus-state.c taus.h taus-gen.c uni-state.c uni.h uni-gen.c uni32-state.c uni32.h uni32-gen.c zuf-state.c zuf.h zuf-gen.c rand-state.c rand.h rand-gen.c 