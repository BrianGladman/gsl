lib_LTLIBRARIES = libgslcblas.la

pkginclude_HEADERS = gsl_cblas.h

INCLUDES= -I$(top_srcdir)

libgslcblas_la_SOURCES = sasum.c saxpy.c scasum.c scnrm2.c scopy.c sdot.c sdsdot.c sgbmv.c sgemm.c sgemv.c sger.c snrm2.c srot.c srotg.c srotm.c srotmg.c ssbmv.c sscal.c sspmv.c sspr.c sspr2.c sswap.c ssymm.c ssymv.c ssyr.c ssyr2.c ssyr2k.c ssyrk.c stbmv.c stbsv.c stpmv.c stpsv.c strmm.c strmv.c strsm.c strsv.c dasum.c daxpy.c dcopy.c ddot.c dgbmv.c dgemm.c dgemv.c dger.c dnrm2.c drot.c drotg.c drotm.c drotmg.c dsbmv.c dscal.c dsdot.c dspmv.c dspr.c dspr2.c dswap.c dsymm.c dsymv.c dsyr.c dsyr2.c dsyr2k.c dsyrk.c dtbmv.c dtbsv.c dtpmv.c dtpsv.c dtrmm.c dtrmv.c dtrsm.c dtrsv.c dzasum.c dznrm2.c caxpy.c ccopy.c cdotc_sub.c cdotu_sub.c cgbmv.c cgemm.c cgemv.c cgerc.c cgeru.c chbmv.c chemm.c chemv.c cher.c cher2.c cher2k.c cherk.c chpmv.c chpr.c chpr2.c cscal.c csscal.c cswap.c csymm.c csyr2k.c csyrk.c ctbmv.c ctbsv.c ctpmv.c ctpsv.c ctrmm.c ctrmv.c ctrsm.c ctrsv.c zaxpy.c zcopy.c zdotc_sub.c zdotu_sub.c zdscal.c zgbmv.c zgemm.c zgemv.c zgerc.c zgeru.c zhbmv.c zhemm.c zhemv.c zher.c zher2.c zher2k.c zherk.c zhpmv.c zhpr.c zhpr2.c zscal.c zswap.c zsymm.c zsyr2k.c zsyrk.c ztbmv.c ztbsv.c ztpmv.c ztpsv.c ztrmm.c ztrmv.c ztrsm.c ztrsv.c

noinst_HEADERS = cblas.h source_asum_c.h source_asum_r.h source_axpy_c.h source_axpy_r.h source_copy_c.h source_copy_r.h source_dot_c.h source_dot_r.h source_gbmv_c.h source_gbmv_r.h source_gemm_c.h source_gemm_r.h source_gemv_c.h source_gemv_r.h source_ger.h source_gerc.h source_geru.h source_hbmv.h source_hemm.h source_hemv.h source_her.h source_her2.h source_her2k.h source_herk.h source_hpmv.h source_hpr.h source_hpr2.h source_iamax_c.h source_iamax_r.h source_nrm2_c.h source_nrm2_r.h source_rot.h source_rotg.h source_rotm.h source_rotmg.h source_sbmv.h source_scal_c.h source_scal_c_s.h source_scal_r.h source_spmv.h source_spr.h source_spr2.h source_swap_c.h source_swap_r.h source_symm_c.h source_symm_r.h source_symv.h source_syr.h source_syr2.h source_syr2k_c.h source_syr2k_r.h source_syrk_c.h source_syrk_r.h source_tbmv_c.h source_tbmv_r.h source_tbsv_c.h source_tbsv_r.h source_tpmv_c.h source_tpmv_r.h source_tpsv_c.h source_tpsv_r.h source_trmm_c.h source_trmm_r.h source_trmv_c.h source_trmv_r.h source_trsm_c.h source_trsm_r.h source_trsv_c.h source_trsv_r.h

check_PROGRAMS = test
TESTS = test

test_LDADD = libgslcblas.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la
test_SOURCES = test.c 
test_INCLUDES = tests.c
