#####################################################################

## I have removed the "automagic" making of Makefile.am from
## Makefile.ami  --- It was just too finicky.  

## If you consider yourself a casual user, then you shouldn't need to
## edit the Makefile.am, but if you do, nobody will mind.

## If instead you think of yourself as a developer/contributor [eg if
## you are working from a distribution you got from the CVS master
## source instead of a numbered release], and if you want to edit the
## Makefile.am file, you should instead edit the Makefile.ami template
## file, and/or modify the 'makeam' shell script.  Then you can remake
## the Makefile.am file using the 'makeam' script:

##      ./makeam Makefile.ami > Makefile.am

## After you are done, be sure the 'cvs commit' the Makefile.am file,
## as well as your changes to 'Makefile.ami' and 'makeam'

## A similar warning applies to the sources *-state.c, *-gen.c, and
## some of the *.h files; if you want to edit them, you should instead
## edit the template files, xxx-state.xc, xxx-gen.xc, and gsl_ran.h
## and then run './makesrc --all' to remake all the built-sources.

#####################################################################

## Some *.h files are in neighboring directories
INCLUDES = -I$(srcdir)/../err -I$(srcdir)/../utils

## Include here the list of algorithms for random number generation
## The first one should be the "favorite" (ie, default)

ALLRAN= taus mrg cmrg uni uni32 zuf rand

include_HEADERS= gsl_ran.h gsl_ran_switch.h

test_switch_SOURCES= test-switch.c switch.c 

switch.c gsl_ran_switch.h:
	$(srcdir)/makesrc $(srcdir) --switch

## The user doesn't really _need_ these, assuming that the
## distribution includes all the .c and .h sources.
## But it is fair for the user to have them, in case one of
## those sources needs to be rebuilt, and also because it
## permits the user to contribute to the project!!
##
EXTRA_DIST= xxx-switch.xc xxx-gen.xc xxx-state.xc \
	makeam makesrc Makefile.ami 

CLEANFILES= tmp.out tmp.in







bin_PROGRAMS= taus test-taus bench-taus mrg test-mrg bench-mrg cmrg test-cmrg bench-cmrg uni test-uni bench-uni uni32 test-uni32 bench-uni32 zuf test-zuf bench-zuf rand test-rand bench-rand test-switch 

lib_LIBRARIES= libgslrandom.a libgslrantaus.a libgslranmrg.a libgslrancmrg.a libgslranuni.a libgslranuni32.a libgslranzuf.a libgslranrand.a 

libgslrantaus_a_SOURCES= taus.c taus-gen.c 
libgslranmrg_a_SOURCES= mrg.c mrg-gen.c 
libgslrancmrg_a_SOURCES= cmrg.c cmrg-gen.c 
libgslranuni_a_SOURCES= uni.c uni-gen.c 
libgslranuni32_a_SOURCES= uni32.c uni32-gen.c 
libgslranzuf_a_SOURCES= zuf.c zuf-gen.c 
libgslranrand_a_SOURCES= rand.c rand-gen.c 

libgslrandom_a_SOURCES= switch.c taus.c mrg.c cmrg.c uni.c uni32.c zuf.c rand.c 

test_switch_DEPENDENCIES= gsl_ran_switch.h taus.h taus.o mrg.h mrg.o cmrg.h cmrg.o uni.h uni.o uni32.h uni32.o zuf.h zuf.o rand.h rand.o 
test_switch_LDADD= taus.o mrg.o cmrg.o uni.o uni32.o zuf.o rand.o 
taus_SOURCES=             ran.c taus.c
taus_DEPENDENCIES=        taus-state.c taus.h taus-gen.o
taus_LDADD=               taus-gen.o
test_taus_SOURCES=        testran.c
test_taus_DEPENDENCIES=   taus.o taus-gen.o taus.h
test_taus_LDADD=          taus.o taus-gen.o ../err/libgsltest.a ../utils/libutils.a
bench_taus_SOURCES=       benchran.c
bench_taus_DEPENDENCIES=  taus.o taus-gen.o taus.h
bench_taus_LDADD=         taus.o taus-gen.o
mrg_SOURCES=             ran.c mrg.c
mrg_DEPENDENCIES=        mrg-state.c mrg.h mrg-gen.o
mrg_LDADD=               mrg-gen.o
test_mrg_SOURCES=        testran.c
test_mrg_DEPENDENCIES=   mrg.o mrg-gen.o mrg.h
test_mrg_LDADD=          mrg.o mrg-gen.o ../err/libgsltest.a ../utils/libutils.a
bench_mrg_SOURCES=       benchran.c
bench_mrg_DEPENDENCIES=  mrg.o mrg-gen.o mrg.h
bench_mrg_LDADD=         mrg.o mrg-gen.o
cmrg_SOURCES=             ran.c cmrg.c
cmrg_DEPENDENCIES=        cmrg-state.c cmrg.h cmrg-gen.o
cmrg_LDADD=               cmrg-gen.o
test_cmrg_SOURCES=        testran.c
test_cmrg_DEPENDENCIES=   cmrg.o cmrg-gen.o cmrg.h
test_cmrg_LDADD=          cmrg.o cmrg-gen.o ../err/libgsltest.a ../utils/libutils.a
bench_cmrg_SOURCES=       benchran.c
bench_cmrg_DEPENDENCIES=  cmrg.o cmrg-gen.o cmrg.h
bench_cmrg_LDADD=         cmrg.o cmrg-gen.o
uni_SOURCES=             ran.c uni.c
uni_DEPENDENCIES=        uni-state.c uni.h uni-gen.o
uni_LDADD=               uni-gen.o
test_uni_SOURCES=        testran.c
test_uni_DEPENDENCIES=   uni.o uni-gen.o uni.h
test_uni_LDADD=          uni.o uni-gen.o ../err/libgsltest.a ../utils/libutils.a
bench_uni_SOURCES=       benchran.c
bench_uni_DEPENDENCIES=  uni.o uni-gen.o uni.h
bench_uni_LDADD=         uni.o uni-gen.o
uni32_SOURCES=             ran.c uni32.c
uni32_DEPENDENCIES=        uni32-state.c uni32.h uni32-gen.o
uni32_LDADD=               uni32-gen.o
test_uni32_SOURCES=        testran.c
test_uni32_DEPENDENCIES=   uni32.o uni32-gen.o uni32.h
test_uni32_LDADD=          uni32.o uni32-gen.o ../err/libgsltest.a ../utils/libutils.a
bench_uni32_SOURCES=       benchran.c
bench_uni32_DEPENDENCIES=  uni32.o uni32-gen.o uni32.h
bench_uni32_LDADD=         uni32.o uni32-gen.o
zuf_SOURCES=             ran.c zuf.c
zuf_DEPENDENCIES=        zuf-state.c zuf.h zuf-gen.o
zuf_LDADD=               zuf-gen.o
test_zuf_SOURCES=        testran.c
test_zuf_DEPENDENCIES=   zuf.o zuf-gen.o zuf.h
test_zuf_LDADD=          zuf.o zuf-gen.o ../err/libgsltest.a ../utils/libutils.a
bench_zuf_SOURCES=       benchran.c
bench_zuf_DEPENDENCIES=  zuf.o zuf-gen.o zuf.h
bench_zuf_LDADD=         zuf.o zuf-gen.o
rand_SOURCES=             ran.c rand.c
rand_DEPENDENCIES=        rand-state.c rand.h rand-gen.o
rand_LDADD=               rand-gen.o
test_rand_SOURCES=        testran.c
test_rand_DEPENDENCIES=   rand.o rand-gen.o rand.h
test_rand_LDADD=          rand.o rand-gen.o ../err/libgsltest.a ../utils/libutils.a
bench_rand_SOURCES=       benchran.c
bench_rand_DEPENDENCIES=  rand.o rand-gen.o rand.h
bench_rand_LDADD=         rand.o rand-gen.o

## Note, BUILT_SOURCES has special meaning to automake
## I'm deliberately trying to avoid having the casual user
## build sources.  If you want to build THE_BUILT_SOURCESS
## then run the program './makesrc'
## THE_BUILT_SOURCES= taus.h taus-state.c taus-gen.c mrg.h mrg-state.c mrg-gen.c cmrg.h cmrg-state.c cmrg-gen.c uni.h uni-state.c uni-gen.c uni32.h uni32-state.c uni32-gen.c zuf.h zuf-state.c zuf-gen.c rand.h rand-state.c rand-gen.c gsl_ran_switch.h switch.c
## If you need to build these sources, use the 'makesrc' script
taus.h taus-state.c taus-gen.c: makesrc
	$(srcdir)/makesrc $(srcdir) taus

mrg.h mrg-state.c mrg-gen.c: makesrc
	$(srcdir)/makesrc $(srcdir) mrg

cmrg.h cmrg-state.c cmrg-gen.c: makesrc
	$(srcdir)/makesrc $(srcdir) cmrg

uni.h uni-state.c uni-gen.c: makesrc
	$(srcdir)/makesrc $(srcdir) uni

uni32.h uni32-state.c uni32-gen.c: makesrc
	$(srcdir)/makesrc $(srcdir) uni32

zuf.h zuf-state.c zuf-gen.c: makesrc
	$(srcdir)/makesrc $(srcdir) zuf

rand.h rand-state.c rand-gen.c: makesrc
	$(srcdir)/makesrc $(srcdir) rand

## This target is kind of like 'distclean'
tooclean:
	/bin/rm -rf .deps
	/bin/rm -f taus.h taus-state.c taus-gen.c
	/bin/rm -f mrg.h mrg-state.c mrg-gen.c
	/bin/rm -f cmrg.h cmrg-state.c cmrg-gen.c
	/bin/rm -f uni.h uni-state.c uni-gen.c
	/bin/rm -f uni32.h uni32-state.c uni32-gen.c
	/bin/rm -f zuf.h zuf-state.c zuf-gen.c
	/bin/rm -f rand.h rand-state.c rand-gen.c
	/bin/rm -f switch.c gsl_ran_switch.h

