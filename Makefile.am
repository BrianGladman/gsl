## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libgsl.a

DEP_FILES=""

SUBDIRS = err ieee-utils vector matrix blas specfunc rng randist fft roots statistics siman integration interpolation histogram monte doc

include_HEADERS = gsl_config.h gsl_complex.h gsl_math.h templates_on.h templates_off.h
noinst_HEADERS = complex_internal.h

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -Werror -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -Dinline='' -g -O4"

# Extend 'make all' through automake to build a single libgsl.a file

# The three following targets work around automake

libgsl.a.o:

libgsl.a.c:

.PHONY: libgsl.a.o libgsl.a.c

# The next target builds a new libgsl.a if any of the subdir/lib.a's
# are newer It has to be kept up to date manually, using the libfixmk.pl
# script.

libgsl.a: err/libgslerr.a err/libgsltest.a ieee-utils/libgslieeeutils.a vector/libgslvector.a matrix/libgslmatrix.a blas/libgslblas.a blas/libgslblasnative.a blas/libgslblascblas.a specfunc/libgslspecfunc.a rng/libgslrng.a randist/libgslrandist.a fft/libgslfft.a roots/libgslroots.a statistics/libgslstatistics.a siman/libgslsiman.a integration/libgslintegration.a interpolation/libgslinterpolation.a histogram/libgslhistogram.a monte/libgslmonte.a
	@rm -f libgsl.a
	$(AR) q libgsl.a err/error.o err/stream.o err/message.o err/strerror.o err/warn.o
	$(AR) q libgsl.a err/test_results.o
	$(AR) q libgsl.a ieee-utils/endian.o ieee-utils/print.o ieee-utils/standardize.o ieee-utils/make_rep.o ieee-utils/env.o ieee-utils/fp.o ieee-utils/read.o
	$(AR) q libgsl.a vector/init.o vector/file.o vector/vector.o
	$(AR) q libgsl.a matrix/init.o matrix/matrix.o matrix/file.o matrix/rowcol.o
	$(AR) q libgsl.a blas/blas.o
	$(AR) q libgsl.a blas/blas_raw_native.o blas/blas_raw_native_cm.o
	$(AR) q libgsl.a blas/blas_raw_cblas.o blas/blas_raw_cblas_cm.o
	$(AR) q libgsl.a specfunc/airy.o specfunc/airy_der.o specfunc/bessel.o specfunc/bessel_I0.o specfunc/bessel_I1.o specfunc/bessel_In.o specfunc/bessel_Inu.o specfunc/bessel_J0.o specfunc/bessel_J1.o specfunc/bessel_Jn.o specfunc/bessel_Jnu.o specfunc/bessel_K0.o specfunc/bessel_K1.o specfunc/bessel_Kn.o specfunc/bessel_Y0.o specfunc/bessel_Y1.o specfunc/bessel_Yn.o specfunc/bessel_Ynu.o specfunc/bessel_amp_phase.o specfunc/bessel_i.o specfunc/bessel_j.o specfunc/bessel_k.o specfunc/bessel_olver.o specfunc/bessel_y.o specfunc/chebyshev.o specfunc/clausen.o specfunc/coupling.o specfunc/coulomb.o specfunc/coulomb_bound.o specfunc/debye.o specfunc/dilog.o specfunc/ellint.o specfunc/elljac.o specfunc/erfc.o specfunc/exp.o specfunc/expint.o specfunc/expint3.o specfunc/fermi_dirac.o specfunc/gegenbauer.o specfunc/gamma.o specfunc/hyperg_0F1.o specfunc/hyperg_2F0.o specfunc/hyperg_1F1.o specfunc/hyperg_2F1.o specfunc/laguerre.o specfunc/legendre_con.o specfunc/legendre_poly.o specfunc/log.o specfunc/poly.o specfunc/pow_int.o specfunc/psi.o specfunc/shint.o specfunc/sig2lamp.o specfunc/sinint.o specfunc/synchrotron.o specfunc/transport.o specfunc/trig.o specfunc/zeta.o
	$(AR) q libgsl.a rng/rng.o rng/default.o rng/cmrg.o rng/slatec.o rng/minstd.o rng/mrg.o rng/mt.o rng/r250.o rng/ran0.o rng/ran1.o rng/ran2.o rng/ran3.o rng/rand.o rng/random.o rng/rand48.o rng/randu.o rng/ranf.o rng/ranlux.o rng/ranmar.o rng/taus.o rng/transputer.o rng/tt.o rng/uni.o rng/uni32.o rng/vax.o rng/zuf.o
	$(AR) q libgsl.a randist/bernoulli.o randist/beta.o randist/binomial.o randist/cauchy.o randist/chisq.o randist/erlang.o randist/exponential.o randist/exppow.o randist/fdist.o randist/flat.o randist/gamma.o randist/gauss.o randist/geometric.o randist/gumbel.o randist/hyperg.o randist/laplace.o randist/levy.o randist/logarithmic.o randist/logistic.o randist/lognormal.o randist/nbinomial.o randist/pareto.o randist/pascal.o randist/poisson.o randist/shuffle.o randist/sphere.o randist/tdist.o randist/weibull.o
	$(AR) q libgsl.a fft/compare.o fft/dft_complex.o fft/c.o fft/c_init.o fft/c_pass_2.o fft/c_pass_3.o fft/c_pass_4.o fft/c_pass_5.o fft/c_pass_6.o fft/c_pass_7.o fft/c_pass_n.o fft/c_radix2.o fft/factorize.o fft/bitreverse.o fft/hc.o fft/hc_init.o fft/hc_pass_2.o fft/hc_pass_3.o fft/hc_pass_4.o fft/hc_pass_5.o fft/hc_pass_n.o fft/hc_radix2.o fft/hc_unpack.o fft/real.o fft/real_init.o fft/real_pass_2.o fft/real_pass_3.o fft/real_pass_4.o fft/real_pass_5.o fft/real_pass_n.o fft/real_radix2.o fft/real_unpack.o
	$(AR) q libgsl.a roots/utility.o roots/bisection.o roots/falsepos.o roots/secant.o roots/newton.o
	$(AR) q libgsl.a statistics/mean.o statistics/variance.o statistics/absdev.o statistics/skew.o statistics/kurtosis.o statistics/lag1.o statistics/p_variance.o statistics/minmax.o statistics/ttest.o statistics/sort.o statistics/median.o statistics/quantiles.o
	$(AR) q libgsl.a siman/siman.o
	$(AR) q libgsl.a integration/integration.o integration/qk15.o integration/qk21.o integration/qk31.o integration/qk41.o integration/qk51.o integration/qk61.o integration/qk.o integration/qng.o integration/qage.o integration/qage_impl.o integration/qpsrt.o integration/qagse.o integration/qagse_impl.o integration/qelg.o integration/err.o integration/workspace.o
	$(AR) q libgsl.a interpolation/accel.o interpolation/akima.o interpolation/bsearch.o interpolation/cspline.o interpolation/interp.o interpolation/linear.o interpolation/tridiag.o
	$(AR) q libgsl.a histogram/add.o histogram/find.o histogram/get.o histogram/init.o histogram/params.o histogram/reset.o histogram/file.o histogram/pdf.o histogram/add2d.o histogram/find2d.o histogram/get2d.o histogram/init2d.o histogram/params2d.o histogram/reset2d.o histogram/file2d.o histogram/pdf2d.o
	$(AR) q libgsl.a monte/vegas.o monte/vegas-print.o monte/miser.o monte/plain.o monte/utils.o monte/init_source.o
	$(RANLIB) libgsl.a
