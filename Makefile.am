## Process this file with automake to produce Makefile.in

# AUTOMAKE_OPTIONS = readme-alpha

SUBDIRS = utils sys err block vector matrix ieee-utils min multimin linalg eigen blas specfunc dht rng randist fft poly roots multiroots statistics siman sum integration interpolation histogram monte ode-initval doc

bin_SCRIPTS = gsl-config

pkginclude_HEADERS = gsl_complex.h gsl_math.h gsl_machine.h gsl_mode.h gsl_precision.h
noinst_HEADERS = complex_internal.h templates_on.h templates_off.h
EXTRA_DIST = autogen.sh gsl-config.in configure.in THANKS KNOWN-PROBLEMS

pkglib_LIBRARIES = libgsl.a

# libgsl_a_LIBADD = blas/libgslblas.a blas/libgslblasnative.a block/libgslblock.a dht/libgsldht.a eigen/libgsleigen.a err/libgslerr.a err/libgsltest.a fft/libgslfft.a histogram/libgslhistogram.a ieee-utils/libgslieeeutils.a integration/libgslintegration.a interpolation/libgslinterpolation.a linalg/libgsllinalg.a matrix/libgslmatrix.a min/libgslmin.a monte/libgslmonte.a multimin/libgslmultimin.a multiroots/libgslmultiroots.a ode-initval/libgslodeiv.a poly/libgslpoly.a randist/libgslrandist.a rng/libgslrng.a roots/libgslroots.a siman/libgslsiman.a specfunc/libgslspecfunc.a statistics/libgslstatistics.a sum/libgslsum.a sys/libgslsys.a utils/libutils.a vector/libgslvector.a

libgsl_a_SOURCES = version.c

