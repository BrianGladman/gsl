noinst_LTLIBRARIES = libgslblas.la
lib_LTLIBRARIES = libgslblasnative.la libgslblascblas.la

check_PROGRAMS = test

pkginclude_HEADERS = 		\
	gsl_blas.h		\
	gsl_blas_types.h	\
	gsl_blas_raw.h		\
	gsl_blas_raw_L1.h	\
	gsl_blas_raw_L2.h	\
	gsl_blas_raw_L3.h	\
	gsl_blas_raw_cm.h	\
	gsl_cblas.h

noinst_HEADERS = matrix_access.h \
	source_asum_c.h		\
	source_asum_r.h		\
	source_axpy_c.h		\
	source_axpy_r.h		\
	source_copy_c.h		\
	source_copy_r.h		\
	source_dot_c.h		\
	source_dot_r.h		\
	source_gbmv_c.h		\
	source_gbmv_r.h		\
	source_gemv_c.h		\
	source_gemv_r.h		\
	source_ger.h		\
	source_gerc.h		\
	source_geru.h		\
	source_hbmv.h		\
	source_her.h		\
	source_her2.h		\
	source_hpmv.h		\
	source_hpr.h		\
	source_hpr2.h		\
	source_iamax_c.h	\
	source_iamax_r.h	\
	source_nrm2_c.h		\
	source_nrm2_r.h		\
	source_rot.h		\
	source_rotg.h		\
	source_rotm.h		\
	source_rotmg.h		\
	source_sbmv.h		\
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
	source_syr2.h		\
	source_tXsv_c.h		\
	source_tXsv_r.h		\
	source_tbmv_c.h		\
	source_tbmv_r.h		\
	source_tbsv_c.h		\
	source_tbsv_r.h		\
	source_tpmv_c.h		\
	source_tpmv_r.h		\
	source_tpsv_c.h         \
	source_tpsv_r.h


INCLUDES= -I$(top_srcdir)

TESTS = test

test_LDADD = libgslblasnative.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la
test_SOURCES = test_blas_raw.c test_cases.c test_cases.h

libgslblas_la_SOURCES = 		\
	blas.c

libgslblasnative_la_SOURCES = 	\
	blas_raw_native.c	\
	blas_raw_native_cm.c

libgslblascblas_la_SOURCES =	\
	blas_raw_cblas.c	\
	blas_raw_cblas_cm.c

