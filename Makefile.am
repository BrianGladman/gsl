## Process this file with automake to produce Makefile.in

# AUTOMAKE_OPTIONS = readme-alpha

SUBDIRS = gsl utils sys test err const complex block vector matrix permutation sort ieee-utils blas linalg eigen specfunc dht qrng rng randist fft poly fit multifit statistics siman sum integration interpolation histogram monte ode-initval roots multiroots min multimin doc

SUBLIBS = block/libgslblock.la complex/libgslcomplex.la dht/libgsldht.la eigen/libgsleigen.la err/libgslerr.la test/libgsltest.la fft/libgslfft.la histogram/libgslhistogram.la ieee-utils/libgslieeeutils.la integration/libgslintegration.la interpolation/libgslinterpolation.la linalg/libgsllinalg.la matrix/libgslmatrix.la min/libgslmin.la monte/libgslmonte.la multimin/libgslmultimin.la multiroots/libgslmultiroots.la ode-initval/libgslodeiv.la permutation/libgslpermutation.la poly/libgslpoly.la fit/libgslfit.la multifit/libgslmultifit.la qrng/libgslqrng.la randist/libgslrandist.la rng/libgslrng.la roots/libgslroots.la siman/libgslsiman.la sort/libgslsort.la specfunc/libgslspecfunc.la statistics/libgslstatistics.la sum/libgslsum.la sys/libgslsys.la utils/libutils.la vector/libgslvector.la

bin_SCRIPTS = gsl-config

pkginclude_HEADERS = gsl_math.h gsl_machine.h gsl_mode.h gsl_precision.h
EXTRA_DIST = autogen.sh gsl-config.in configure.in THANKS KNOWN-PROBLEMS gsl.spec gsl.m4 MACHINES

lib_LTLIBRARIES = libgsl.la
libgsl_la_SOURCES = version.c
libgsl_la_LIBADD = $(SUBLIBS)
noinst_HEADERS = complex_internal.h templates_on.h templates_off.h

m4datadir = $(datadir)/aclocal
m4data_DATA = gsl.m4

#bin_PROGRAMS = main dummy
#dummy_SOURCES = version.c
#dummy_LDADD = $(SUBLIBS)
#main_SOURCES = version.c
#main_LDADD = libgsl.la
