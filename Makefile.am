## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libgsl.a

DEP_FILES = ""

SUBDIRS = err sys utils ieee-utils vector matrix linalg blas specfunc rng randist fft poly roots statistics siman sum integration interpolation histogram monte doc

include_HEADERS = gsl_complex.h gsl_math.h gsl_machine.h gsl_mode.h gsl_precision.h
noinst_HEADERS = complex_internal.h templates_on.h templates_off.h
EXTRA_DIST = autogen.sh

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

libgsl.a: err/libgslerr.a err/libgsltest.a sys/libgslsys.a ieee-utils/libgslieeeutils.a vector/libgslvector.a matrix/libgslmatrix.a linalg/libgsllinalg.a blas/libgslblas.a blas/libgslblasnative.a blas/libgslblascblas.a specfunc/libgslspecfunc.a rng/libgslrng.a randist/libgslrandist.a fft/libgslfft.a poly/libgslpoly.a roots/libgslroots.a statistics/libgslstatistics.a siman/libgslsiman.a sum/libgslsum.a integration/libgslintegration.a interpolation/libgslinterpolation.a histogram/libgslhistogram.a monte/libgslmonte.a
	@rm -f libgsl.a
	$(AR) q libgsl.a err/error.o err/stream.o err/message.o err/strerror.o err/warn.o
	$(AR) q libgsl.a err/test_results.o
	$(AR) q libgsl.a sys/minmax.o sys/prec.o
	$(AR) q libgsl.a ieee-utils/endian.o ieee-utils/print.o ieee-utils/standardize.o ieee-utils/make_rep.o ieee-utils/env.o ieee-utils/fp.o ieee-utils/read.o
	$(AR) q libgsl.a vector/init.o vector/file.o vector/vector.o
	$(AR) q libgsl.a matrix/init.o matrix/matrix.o matrix/file.o matrix/rowcol.o
	$(AR) q libgsl.a linalg/eigen_sort.o linalg/linalg_simple.o linalg/multiply.o
	$(AR) q libgsl.a blas/blas.o
	$(AR) q libgsl.a blas/blas_raw_native.o blas/blas_raw_native_cm.o
	$(AR) q libgsl.a blas/blas_raw_cblas.o blas/blas_raw_cblas_cm.o
	$(AR) q libgsl.a specfunc/airy.o specfunc/airy_der.o specfunc/atanint.o specfunc/bessel.o specfunc/bessel_I0.o specfunc/bessel_I1.o specfunc/bessel_In.o specfunc/bessel_Inu.o specfunc/bessel_J0.o specfunc/bessel_J1.o specfunc/bessel_Jn.o specfunc/bessel_Jnu.o specfunc/bessel_K0.o specfunc/bessel_K1.o specfunc/bessel_Kn.o specfunc/bessel_Knu.o specfunc/bessel_Y0.o specfunc/bessel_Y1.o specfunc/bessel_Yn.o specfunc/bessel_Ynu.o specfunc/bessel_amp_phase.o specfunc/bessel_i.o specfunc/bessel_j.o specfunc/bessel_k.o specfunc/bessel_olver.o specfunc/bessel_temme.o specfunc/bessel_y.o specfunc/beta.o specfunc/chebyshev.o specfunc/clausen.o specfunc/coulomb.o specfunc/coupling.o specfunc/coulomb_bound.o specfunc/dawson.o specfunc/debye.o specfunc/dilog.o specfunc/elementary.o specfunc/ellint.o specfunc/elljac.o specfunc/erfc.o specfunc/exp.o specfunc/expint.o specfunc/expint3.o specfunc/fermi_dirac.o specfunc/gegenbauer.o specfunc/gamma.o specfunc/gamma_inc.o specfunc/hyperg_0F1.o specfunc/hyperg_2F0.o specfunc/hyperg_1F1.o specfunc/hyperg_2F1.o specfunc/hyperg_U.o specfunc/hyperg.o specfunc/laguerre.o specfunc/legendre_H3d.o specfunc/legendre_Qn.o specfunc/legendre_con.o specfunc/legendre_poly.o specfunc/log.o specfunc/poch.o specfunc/poly.o specfunc/pow_int.o specfunc/psi.o specfunc/shint.o specfunc/sinint.o specfunc/synchrotron.o specfunc/transport.o specfunc/trig.o specfunc/zeta.o
	$(AR) q libgsl.a rng/rng.o rng/default.o rng/cmrg.o rng/gfsr4.o rng/slatec.o rng/minstd.o rng/mrg.o rng/mt.o rng/r250.o rng/ran0.o rng/ran1.o rng/ran2.o rng/ran3.o rng/rand.o rng/random.o rng/rand48.o rng/randu.o rng/ranf.o rng/ranlux.o rng/ranmar.o rng/taus.o rng/transputer.o rng/tt.o rng/uni.o rng/uni32.o rng/vax.o rng/zuf.o
	$(AR) q libgsl.a randist/bernoulli.o randist/beta.o randist/binomial.o randist/cauchy.o randist/chisq.o randist/discrete.o randist/erlang.o randist/exponential.o randist/exppow.o randist/fdist.o randist/flat.o randist/gamma.o randist/gauss.o randist/geometric.o randist/gumbel.o randist/hyperg.o randist/laplace.o randist/levy.o randist/logarithmic.o randist/logistic.o randist/lognormal.o randist/nbinomial.o randist/pareto.o randist/pascal.o randist/poisson.o randist/rayleigh.o randist/shuffle.o randist/sphere.o randist/tdist.o randist/weibull.o
	$(AR) q libgsl.a fft/compare.o fft/dft.o fft/c.o fft/c_float.o fft/factorize.o fft/hc.o fft/real.o
	$(AR) q libgsl.a poly/solve_quadratic.o poly/solve_cubic.o poly/zsolve_quadratic.o poly/zsolve_cubic.o
	$(AR) q libgsl.a roots/bisection.o roots/brent.o roots/falsepos.o roots/newton.o roots/secant.o roots/steffenson.o roots/convergence.o roots/fsolver.o roots/fdfsolver.o
	$(AR) q libgsl.a statistics/mean.o statistics/variance.o statistics/absdev.o statistics/skew.o statistics/kurtosis.o statistics/lag1.o statistics/p_variance.o statistics/minmax.o statistics/ttest.o statistics/sort.o statistics/median.o statistics/quantiles.o
	$(AR) q libgsl.a siman/siman.o
	$(AR) q libgsl.a sum/levin_u.o sum/levin_uerr.o
	$(AR) q libgsl.a integration/integration.o integration/qk15.o integration/qk21.o integration/qk31.o integration/qk41.o integration/qk51.o integration/qk61.o integration/qk.o integration/qng.o integration/qage.o integration/qage_impl.o integration/qpsrt.o integration/qagse.o integration/qagse_impl.o integration/qelg.o integration/err.o integration/workspace.o
	$(AR) q libgsl.a interpolation/accel.o interpolation/akima.o interpolation/bsearch.o interpolation/cspline.o interpolation/interp.o interpolation/linear.o interpolation/tridiag.o
	$(AR) q libgsl.a histogram/add.o histogram/find.o histogram/get.o histogram/init.o histogram/params.o histogram/reset.o histogram/file.o histogram/pdf.o histogram/add2d.o histogram/find2d.o histogram/get2d.o histogram/init2d.o histogram/params2d.o histogram/reset2d.o histogram/file2d.o histogram/pdf2d.o
	$(AR) q libgsl.a monte/vegas.o monte/vegas-print.o monte/miser.o monte/plain.o monte/utils.o
	$(RANLIB) libgsl.a
