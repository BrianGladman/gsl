pkglib_LIBRARIES = 			\
	libgslblas.a 			\
	libgslblasnative.a		\
	libgslblascblas.a

# check_PROGRAMS = test_blas

pkginclude_HEADERS = 		\
	gsl_blas.h		\
	gsl_blas_types.h	\
	gsl_blas_raw.h		\
	gsl_blas_raw_L1.h	\
	gsl_blas_raw_cm.h

noinst_HEADERS =  \
	source_asum_r.h		\
	source_asum_c.h		\
	source_axpy_c.h		\
	source_axpy_r.h		\
	source_copy_c.h		\
	source_copy_r.h		\
	source_dot_c.h		\
	source_dot_r.h		\
	source_gemv.h		\
	source_ger.h		\
	source_iamax_r.h	\
	source_iamax_c.h	\
	source_nrm2_c.h		\
	source_nrm2_r.h		\
	source_rot.h		\
	source_rotg.h		\
	source_scal_c.h		\
	source_scal_c_s.h	\
	source_scal_r.h		\
	source_spmv.h		\
	source_spr.h		\
	source_spr2.h		\
	source_swap_c.h		\
	source_swap_r.h		\
	source_symv.h		\
	source_syr.h		\
	source_syr2.h


INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err -I$(srcdir)/../matrix -I$(srcdir)/../vector -I$(srcdir)/../block

# TESTS = test_blas

#test_blas_LDADD = libgslblas.a ../err/libgslerr.a ../err/libgsltest.a
#test_blas_DEPENDENCIES = libgslblas.a
#test_blas_SOURCES = test_blas.c

libgslblas_a_SOURCES = 		\
	blas.c

libgslblasnative_a_SOURCES = 	\
	blas_raw_native.c	\
	blas_raw_native_cm.c

libgslblascblas_a_SOURCES =	\
	blas_raw_cblas.c	\
	blas_raw_cblas_cm.c


strict:
	$(MAKE) CFLAGS="-ansi -pedantic -W -Wall -Wtraditional -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Waggregate-return -Wstrict-prototypes -fshort-enums -fno-common -Wmissing-prototypes -Wnested-externs -g -O4 -Dinline=''"
