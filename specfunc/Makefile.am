lib_LIBRARIES = libgslspecfunc.a 

check_PROGRAMS = test_sf

include_HEADERS = \
	gsl_sf.h 	    	\
	gsl_specfunc.h		\
	gsl_sf_airy.h       	\
	gsl_sf_bessel.h     	\
	gsl_sf_chebyshev.h  	\
	gsl_sf_coulomb.h	\
	gsl_sf_dilog.h      	\
	gsl_sf_erf.h        	\
	gsl_sf_expint.h        	\
	gsl_sf_fermi_dirac.h	\
	gsl_sf_gamma.h      	\
	gsl_sf_laguerre.h   	\
	gsl_sf_legendre.h	\
	gsl_sf_log.h 		\
	gsl_sf_pow_int.h 	\
	gsl_sf_psi.h		\
	gsl_sf_trig.h


noinst_HEADERS = \
	bessel_amp_phase.h	\
	bessel.h		\
	airy_impl.h		\
	bessel_I0_impl.h	\
	bessel_In_impl.h	\
	bessel_J0_impl.h	\
	bessel_J1_impl.h	\
	bessel_Jn_impl.h	\
	bessel_K0_impl.h	\
	bessel_K1_impl.h	\
	bessel_Kn_impl.h	\
	bessel_Y0_impl.h	\
	bessel_Y1_impl.h	\
	expint_impl.h		\
	shint_impl.h		\
	sinint_impl.h
	

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err -I$(srcdir)/../integration 

TESTS = test_sf

test_sf_LDADD = libgslspecfunc.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_sf_DEPENDENCIES = libgslspecfunc.a
test_sf_SOURCES = test_sf.c

libgslspecfunc_a_SOURCES = 	\
	airy.c			\
	airy_der.c		\
	bessel.c		\
	bessel.h                \
	bessel_I0.c		\
	bessel_I0_impl.h	\
	bessel_I1.c		\
	bessel_I1_impl.h	\
	bessel_In.c		\
	bessel_In_impl.h	\
	bessel_J0.c		\
	bessel_J0_impl.h	\
	bessel_J1.c		\
	bessel_J1_impl.h	\
	bessel_Jn.c		\
	bessel_Jn_impl.h	\
	bessel_K0.c		\
	bessel_K0_impl.h	\
	bessel_K1.c		\
	bessel_K1_impl.h	\
	bessel_Kn.c		\
	bessel_Kn_impl.h	\
	bessel_Y0.c		\
	bessel_Y0_impl.h	\
	bessel_Y1.c		\
	bessel_Y1_impl.h	\
	bessel_Yn.c		\
	bessel_amp_phase.c	\
        bessel_amp_phase.h      \
	bessel_j.c		\
	chebyshev.c		\
#	coulomb.c		\ doesn't work yet
	dilog.c			\
	erfc.c			\
	expint.c		\
	fermi_dirac.c		\
	gamma.c			\
	laguerre.c		\
	legendre.c		\
	legendre_assoc.c	\
	log.c			\
	pow_int.c		\
	psi.c			\
        recurse.h               \
	shint.c			\
	expint_impl.h		\
	shint_impl.h		\
	sig2lamp.c		\
	sinint.c		\
	sinint_impl.h		\
	t-erf.c			\
	trig.c

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4"
