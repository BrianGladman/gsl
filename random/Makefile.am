
#####################################################################

## I have removed the "automagic" making of Makefile.lam from
## Makefile.lami  --- It was just too finicky.  

## If you consider yourself a casual user, then you shouldn't need to
## edit the Makefile.lam, but if you do, nobody will mind.

## If instead you think of yourself as a developer/contributor [eg if
## you are working from a distribution you got from the CVS master
## source instead of a numbered release], and if you want to edit the
## Makefile.lam file, you should instead edit the Makefile.lami template
## file, and/or modify the 'makeam' shell script.  Then you can remake
## the Makefile.lam file using the 'makeam' script:

##      ./makeam Makefile.lami > Makefile.lam

## After you are done, be sure the 'cvs commit' the Makefile.lam file,
## as well as your changes to 'Makefile.lami' and 'makeam'

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
test_switch_LDADD=../err/libgsltest.la ../utils/libutils.a

switch.c gsl_ran_switch.h:
	$(srcdir)/makesrc --switch

## The user doesn't really _need_ these, assuming that the
## distribution includes all the .c and .h sources.
## But it is fair for the user to have them, in case one of
## those sources needs to be rebuilt, and also because it
## permits the user to contribute to the project!!
##
EXTRA_DIST= xxx-switch.xc xxx-gen.xc xxx-state.xc \
	makeam makesrc Makefile.lami 

CLEANFILES= tmp.lout tmp.in







bin_PROGRAMS= taus mrg cmrg uni uni32 zuf rand test-switch


check_PROGRAMS= test-taus bench-taus test-mrg bench-mrg test-cmrg bench-cmrg test-uni bench-uni test-uni32 bench-uni32 test-zuf bench-zuf test-rand bench-rand


lib_LTLIBRARIES= libgslrandom.la  libgslrantaus.la libgslranmrg.la libgslrancmrg.la libgslranuni.la libgslranuni32.la libgslranzuf.la libgslranrand.la


libgslrantaus_la_SOURCES= taus.c taus-gen.c 
libgslranmrg_la_SOURCES= mrg.c mrg-gen.c 
libgslrancmrg_la_SOURCES= cmrg.c cmrg-gen.c 
libgslranuni_la_SOURCES= uni.c uni-gen.c 
libgslranuni32_la_SOURCES= uni32.c uni32-gen.c 
libgslranzuf_la_SOURCES= zuf.c zuf-gen.c 
libgslranrand_la_SOURCES= rand.c rand-gen.c 

libgslrandom_la_SOURCES=switch.c  taus.c mrg.c cmrg.c uni.c uni32.c zuf.c rand.c


TESTS= test-taus test-mrg test-cmrg test-uni test-uni32 test-zuf test-rand


test_switch_DEPENDENCIES=gsl_ran_switch.h taus.h taus.lo mrg.h mrg.lo cmrg.h cmrg.lo uni.h uni.lo uni32.h uni32.lo zuf.h zuf.lo rand.h rand.lo

test_switch_LDADD= taus.lo mrg.lo cmrg.lo uni.lo uni32.lo zuf.lo rand.lo ../utils/libutils.a 

taus_SOURCES=             ran.c taus.c taus.h
EXTRA_taus_SOURCES=       taus-state.c
taus_DEPENDENCIES=        taus-state.c taus-gen.lo
taus_LDADD=               taus-gen.lo
test_taus_SOURCES=        testran.c
test_taus_DEPENDENCIES=   taus.lo taus-gen.lo taus.h
test_taus_LDADD=          taus.lo taus-gen.lo ../err/libgsltest.la ../utils/libutils.a
bench_taus_SOURCES=       benchran.c
bench_taus_DEPENDENCIES=  taus.lo taus-gen.lo taus.h
bench_taus_LDADD=         taus.lo taus-gen.lo
mrg_SOURCES=             ran.c mrg.c mrg.h
EXTRA_mrg_SOURCES=       mrg-state.c
mrg_DEPENDENCIES=        mrg-state.c mrg-gen.lo
mrg_LDADD=               mrg-gen.lo
test_mrg_SOURCES=        testran.c
test_mrg_DEPENDENCIES=   mrg.lo mrg-gen.lo mrg.h
test_mrg_LDADD=          mrg.lo mrg-gen.lo ../err/libgsltest.la ../utils/libutils.a
bench_mrg_SOURCES=       benchran.c
bench_mrg_DEPENDENCIES=  mrg.lo mrg-gen.lo mrg.h
bench_mrg_LDADD=         mrg.lo mrg-gen.lo
cmrg_SOURCES=             ran.c cmrg.c cmrg.h
EXTRA_cmrg_SOURCES=       cmrg-state.c
cmrg_DEPENDENCIES=        cmrg-state.c cmrg-gen.lo
cmrg_LDADD=               cmrg-gen.lo
test_cmrg_SOURCES=        testran.c
test_cmrg_DEPENDENCIES=   cmrg.lo cmrg-gen.lo cmrg.h
test_cmrg_LDADD=          cmrg.lo cmrg-gen.lo ../err/libgsltest.la ../utils/libutils.a
bench_cmrg_SOURCES=       benchran.c
bench_cmrg_DEPENDENCIES=  cmrg.lo cmrg-gen.lo cmrg.h
bench_cmrg_LDADD=         cmrg.lo cmrg-gen.lo
uni_SOURCES=             ran.c uni.c uni.h
EXTRA_uni_SOURCES=       uni-state.c
uni_DEPENDENCIES=        uni-state.c uni-gen.lo
uni_LDADD=               uni-gen.lo
test_uni_SOURCES=        testran.c
test_uni_DEPENDENCIES=   uni.lo uni-gen.lo uni.h
test_uni_LDADD=          uni.lo uni-gen.lo ../err/libgsltest.la ../utils/libutils.a
bench_uni_SOURCES=       benchran.c
bench_uni_DEPENDENCIES=  uni.lo uni-gen.lo uni.h
bench_uni_LDADD=         uni.lo uni-gen.lo
uni32_SOURCES=             ran.c uni32.c uni32.h
EXTRA_uni32_SOURCES=       uni32-state.c
uni32_DEPENDENCIES=        uni32-state.c uni32-gen.lo
uni32_LDADD=               uni32-gen.lo
test_uni32_SOURCES=        testran.c
test_uni32_DEPENDENCIES=   uni32.lo uni32-gen.lo uni32.h
test_uni32_LDADD=          uni32.lo uni32-gen.lo ../err/libgsltest.la ../utils/libutils.a
bench_uni32_SOURCES=       benchran.c
bench_uni32_DEPENDENCIES=  uni32.lo uni32-gen.lo uni32.h
bench_uni32_LDADD=         uni32.lo uni32-gen.lo
zuf_SOURCES=             ran.c zuf.c zuf.h
EXTRA_zuf_SOURCES=       zuf-state.c
zuf_DEPENDENCIES=        zuf-state.c zuf-gen.lo
zuf_LDADD=               zuf-gen.lo
test_zuf_SOURCES=        testran.c
test_zuf_DEPENDENCIES=   zuf.lo zuf-gen.lo zuf.h
test_zuf_LDADD=          zuf.lo zuf-gen.lo ../err/libgsltest.la ../utils/libutils.a
bench_zuf_SOURCES=       benchran.c
bench_zuf_DEPENDENCIES=  zuf.lo zuf-gen.lo zuf.h
bench_zuf_LDADD=         zuf.lo zuf-gen.lo
rand_SOURCES=             ran.c rand.c rand.h
EXTRA_rand_SOURCES=       rand-state.c
rand_DEPENDENCIES=        rand-state.c rand-gen.lo
rand_LDADD=               rand-gen.lo
test_rand_SOURCES=        testran.c
test_rand_DEPENDENCIES=   rand.lo rand-gen.lo rand.h
test_rand_LDADD=          rand.lo rand-gen.lo ../err/libgsltest.la ../utils/libutils.a
bench_rand_SOURCES=       benchran.c
bench_rand_DEPENDENCIES=  rand.lo rand-gen.lo rand.h
bench_rand_LDADD=         rand.lo rand-gen.lo

## Note, BUILT_SOURCES has special meaning to automake
## I'm deliberately trying to avoid having the casual user
## build sources.  If you want to build THE_BUILT_SOURCESS
## then run the program './makesrc'
## THE_BUILT_SOURCES= taus.h taus-state.c taus-gen.c mrg.h mrg-state.c mrg-gen.c cmrg.h cmrg-state.c cmrg-gen.c uni.h uni-state.c uni-gen.c uni32.h uni32-state.c uni32-gen.c zuf.h zuf-state.c zuf-gen.c rand.h rand-state.c rand-gen.c gsl_ran_switch.h switch.c

## If you need to build these sources, use the 'makesrc' script
# These are commented out now that the files are under CVS 
# The relative dates seem to cause a problem for make distcheck
# taus.h taus-state.c taus-gen.c: makesrc
#	$(srcdir)/makesrc taus

# mrg.h mrg-state.c mrg-gen.c: makesrc
#	$(srcdir)/makesrc mrg

# cmrg.h cmrg-state.c cmrg-gen.c: makesrc
#	$(srcdir)/makesrc cmrg

# uni.h uni-state.c uni-gen.c: makesrc
#	$(srcdir)/makesrc uni

# uni32.h uni32-state.c uni32-gen.c: makesrc
#	$(srcdir)/makesrc uni32

# zuf.h zuf-state.c zuf-gen.c: makesrc
#	$(srcdir)/makesrc zuf

# rand.h rand-state.c rand-gen.c: makesrc
#	$(srcdir)/makesrc rand

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

