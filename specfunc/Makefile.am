lib_LIBRARIES = libgslspecfunc.a 

# noinst_PROGRAMS = 

include_HEADERS = gsl_sf.h gsl_sf_chebyshev.h gsl_sf_gamma.h gsl_sf_pow_int.h \
	gsl_sf_airy.h       gsl_sf_dilog.h      gsl_sf_laguerre.h   gsl_sf_trig.h \
	gsl_sf_bessel.h     gsl_sf_erf.h        gsl_sf_log.h


noinst_HEADERS = bessel_amp_phase.h

INCLUDES = -I$(srcdir)/../err -I$(srcdir)/../integration

TESTS = test_trap test_radix2 test

#test_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
#test_trap_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
#test_radix2_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
#errs_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
#benchmark_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a

#test_DEPENDENCIES = libgslfft.a
#test_trap_DEPENDENCIES = libgslfft.a
#test_radix2_DEPENDENCIES = libgslfft.a
#benchmark_DEPENDENCIES = libgslfft.a

#test_SOURCES = test.c
#test_trap_SOURCES = test_trap.c
#test_radix2_SOURCES = test_radix2.c
#benchmark_SOURCES = benchmark.c 

libgslspecfunc_a_SOURCES = 	\
	airy.c			\
	bessel.c		\
	bessel_I0.c		\
	bessel_I1.c		\
	bessel_J0.c		\
	bessel_J1.c		\
	bessel_K0.c		\
	bessel_K1.c		\
	bessel_Y0.c		\
	bessel_Y1.c		\
	bessel_amp_phase.c	\
	chebyshev.c		\
	dilog.c			\
	erfc.c			\
	gamma.c			\
	laguerre.c		\
	log.c			\
	pow_int.c		\
	sig2lamp.c		\
	t-erf.c			\
	trig.c
