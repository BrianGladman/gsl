noinst_LTLIBRARIES = libgslmultifit.la 

pkginclude_HEADERS = gsl_multifit.h gsl_multifit_nlin.h

INCLUDES= -I$(top_builddir)

libgslmultifit_la_SOURCES = multilinear.c work.c lmder.c fsolver.c fdfsolver.c convergence.c gradient.c covar.c

noinst_HEADERS = lmutil.c lmpar.c lmset.c lmiterate.c qrsolv.c test_filip.c test_fn.c test_longley.c test_pontius.c test_fn.c

check_PROGRAMS = test

TESTS = test

test_SOURCES = test.c
test_LDADD = libgslmultifit.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../blas/libgslblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la
test_DEPENDENCIES = libgslmultifit.la

#noinst_PROGRAMS = demo
#demo_SOURCES = demo.c
#demo_LDADD = libgslmultifit.la ../linalg/libgsllinalg.la ../permutation/libgslpermutation.la ../blas/libgslblas.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../randist/libgslrandist.la ../rng/libgslrng.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la

