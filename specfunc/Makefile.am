lib_LIBRARIES = libgslspecfunc.a 

noinst_PROGRAMS = test_sf

include_HEADERS = \
	gsl_sf.h 	    	\
	gsl_specfunc.h		\
	gsl_sf_airy.h       	\
	gsl_sf_bessel.h     	\
	gsl_sf_chebyshev.h  	\
	gsl_sf_coulomb.h	\
	gsl_sf_dilog.h      	\
	gsl_sf_erf.h        	\
	gsl_sf_fermi_dirac.h	\
	gsl_sf_gamma.h      	\
	gsl_sf_laguerre.h   	\
	gsl_sf_legendre.h	\
	gsl_sf_log.h 		\
	gsl_sf_pow_int.h 	\
	gsl_sf_psi.h		\
	gsl_sf_trig.h


noinst_HEADERS = bessel_amp_phase.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err -I$(srcdir)/../integration 

TESTS = test_sf

test_sf_LDADD = libgslspecfunc.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_sf_DEPENDENCIES = libgslspecfunc.a
test_sf_SOURCES = test_sf.c

libgslspecfunc_a_SOURCES = 	\
	airy.c			\
	bessel.c		\
	bessel_I0.c		\
	bessel_I1.c		\
	bessel_In.c		\
	bessel_J0.c		\
	bessel_J1.c		\
	bessel_Jn.c		\
	bessel_K0.c		\
	bessel_K1.c		\
	bessel_Kn.c		\
	bessel_Y0.c		\
	bessel_Y1.c		\
	bessel_Yn.c		\
	bessel_amp_phase.c	\
	bessel_j.c		\
	chebyshev.c		\
	coulomb.c		\
	dilog.c			\
	erfc.c			\
	fermi_dirac.c		\
	gamma.c			\
	laguerre.c		\
	legendre.c		\
	log.c			\
	pow_int.c		\
	psi.c			\
	sig2lamp.c		\
	t-erf.c			\
	trig.c

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4"
