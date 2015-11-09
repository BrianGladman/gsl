## Process this file with automake to produce Makefile.in

check_PROGRAMS = blas block cblas cdf cheb combination multiset const diff eigen fft fftmr fftreal fitting fitting2 fitting3 fitreg histogram histogram2d ieee ieeeround integration interp interp2d intro linalglu largefit matrix matrixw min monte ntupler ntuplew ode-initval permseq permshuffle polyroots qrng randpoisson randwalk rng rngunif robfit rootnewt roots siman sortsmall specfun specfun_e rstat rquantile stat statsort sum vector vectorr vectorview vectorw dwt nlfit interpp eigen_nonsymm bspline poisson interp_compare spmatrix

examples_src = blas.c block.c cblas.c cdf.c cheb.c combination.c multiset.c const.c demo_fn.c diff.c eigen.c fft.c fftmr.c fftreal.c fitting.c fitting2.c fitting3.c fitreg.c histogram.c histogram2d.c ieee.c ieeeround.c integration.c interp.c interp2d.c intro.c linalglu.c largefit.c matrix.c matrixw.c min.c monte.c ntupler.c ntuplew.c ode-initval.c odefixed.c permseq.c permshuffle.c polyroots.c qrng.c randpoisson.c randwalk.c rng.c rngunif.c robfit.c rootnewt.c roots.c siman.c sortsmall.c specfun.c specfun_e.c rstat.c rquantile.c stat.c statsort.c sum.c vector.c vectorr.c vectorview.c vectorw.c demo_fn.h dwt.c expfit.c nlfit.c interpp.c eigen_nonsymm.c bspline.c multimin.c multiminfn.c nmsimplex.c ode-initval-low-level.c poisson.c interp_compare.c spmatrix.c

examples_txt = blas.txt block.txt cblas.txt cdf.txt combination.txt multiset.txt const.txt diff.txt integration.txt intro.txt largefit.txt largefit2.txt linalglu.txt min.txt polyroots.txt randpoisson.2.txt randpoisson.txt rng.txt rngunif.2.txt rngunif.txt sortsmall.txt specfun.txt specfun_e.txt rquantile.txt stat.txt statsort.txt sum.txt vectorview.txt ecg.dat dwt.dat multimin.txt nmsimplex.txt fitreg.txt interp_compare.txt interp2d.txt

dist_noinst_DATA = $(examples_src) $(examples_txt)

LDADD = ../../libgsl.la ../../cblas/libgslcblas.la
AM_DEFAULT_SOURCE_EXT = .c
