pkglib_LIBRARIES = libgslspecfunc.a 

check_PROGRAMS = test_sf

pkginclude_HEADERS = \
	gsl_sf.h 	    	\
	gsl_specfunc.h		\
	gsl_sf_airy.h       	\
	gsl_sf_bessel.h     	\
	gsl_sf_chebyshev.h  	\
	gsl_sf_clausen.h	\
	gsl_sf_coupling.h	\
	gsl_sf_coulomb.h	\
	gsl_sf_dawson.h		\
	gsl_sf_debye.h      	\
	gsl_sf_dilog.h      	\
	gsl_sf_elementary.h	\
	gsl_sf_ellint.h        	\
	gsl_sf_elljac.h        	\
	gsl_sf_erf.h        	\
	gsl_sf_exp.h        	\
	gsl_sf_expint.h        	\
	gsl_sf_fermi_dirac.h	\
	gsl_sf_gegenbauer.h	\
	gsl_sf_gamma.h      	\
	gsl_sf_hyperg.h		\
	gsl_sf_laguerre.h   	\
	gsl_sf_legendre.h	\
	gsl_sf_log.h 		\
	gsl_sf_poly.h		\
	gsl_sf_pow_int.h 	\
	gsl_sf_psi.h		\
        gsl_sf_result.h         \
	gsl_sf_synchrotron.h	\
	gsl_sf_transport.h	\
	gsl_sf_trig.h		\
	gsl_sf_zeta.h


noinst_HEADERS = \
	bessel_amp_phase.h	\
	bessel_olver.h		\
	bessel_temme.h		\
	bessel.h		\
	hyperg.h		\
	legendre.h


INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err 

TESTS = test_sf

test_sf_LDADD = libgslspecfunc.a ../sys/libgslsys.a \
	../err/libgslerr.a ../err/libgsltest.a
test_sf_DEPENDENCIES = libgslspecfunc.a
test_sf_SOURCES = test_sf.c		\
		  test_sf.h		\
		  test_airy.c		\
		  test_bessel.c		\
		  test_dilog.c		\
		  test_gamma.c		\
		  test_hyperg.c		\
		  test_legendre.c


libgslspecfunc_a_SOURCES = 	\
	airy.c			\
	airy_der.c		\
	airy_zero.c		\
	atanint.c		\
	bessel.c		\
	bessel.h                \
	bessel_I0.c		\
	bessel_I1.c		\
	bessel_In.c		\
	bessel_Inu.c		\
	bessel_J0.c		\
	bessel_J1.c		\
	bessel_Jn.c		\
	bessel_Jnu.c		\
	bessel_K0.c		\
	bessel_K1.c		\
	bessel_Kn.c		\
	bessel_Knu.c		\
	bessel_Y0.c		\
	bessel_Y1.c		\
	bessel_Yn.c		\
	bessel_Ynu.c		\
	bessel_amp_phase.c	\
        bessel_amp_phase.h      \
	bessel_i.c		\
	bessel_j.c		\
	bessel_k.c		\
	bessel_olver.c		\
	bessel_temme.c		\
	bessel_y.c		\
	bessel_zero.c		\
	bessel_sequence.c	\
	beta.c			\
	chebyshev.c		\
	clausen.c		\
	coulomb.c		\
	coupling.c		\
	coulomb_bound.c		\
	dawson.c		\
	debye.c			\
	dilog.c			\
	elementary.c		\
	ellint.c		\
	elljac.c		\
	erfc.c			\
	exp.c			\
	expint.c		\
	expint3.c		\
	fermi_dirac.c		\
	gegenbauer.c		\
	gamma.c			\
	gamma_inc.c		\
	hyperg_0F1.c		\
	hyperg_2F0.c		\
	hyperg_1F1.c		\
	hyperg_2F1.c		\
	hyperg_U.c		\
	hyperg.c		\
	laguerre.c		\
	legendre_H3d.c		\
	legendre_Qn.c		\
	legendre_con.c		\
	legendre_poly.c		\
	log.c			\
	poch.c			\
	poly.c			\
	pow_int.c		\
	psi.c			\
        recurse.h               \
	shint.c			\
	sinint.c		\
	synchrotron.c		\
	transport.c		\
	trig.c			\
	zeta.c

strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
