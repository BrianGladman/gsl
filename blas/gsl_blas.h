/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_BLAS_H_
#define GSL_BLAS_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


/* ========================================================================
 * Level 1
 * ========================================================================
 */

int gsl_blas_sdsdot (float alpha,
                     const gsl_vector_float * X,
                     const gsl_vector_float * Y,
                     float * result
                     );

int gsl_blas_dsdot (const gsl_vector_float * X,
                    const gsl_vector_float * Y,
                    double * result
                    );

int gsl_blas_sdot (const gsl_vector_float * X,
                   const gsl_vector_float * Y,
		   float * result
		   );

int gsl_blas_ddot (const gsl_vector * X,
                   const gsl_vector * Y,
		   double * result
		   );


int  gsl_blas_cdotu (const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_complex_float * dotu);

int  gsl_blas_cdotc (const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_complex_float * dotc);

int  gsl_blas_zdotu (const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_complex * dotu);

int  gsl_blas_zdotc (const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_complex * dotc);


float  gsl_blas_snrm2  (const gsl_vector_float * X);
float  gsl_blas_sasum  (const gsl_vector_float * X);
double gsl_blas_dnrm2  (const gsl_vector * X);
double gsl_blas_dasum  (const gsl_vector * X);
float  gsl_blas_scnrm2 (const gsl_vector_complex_float * X);
float  gsl_blas_scasum (const gsl_vector_complex_float * X);
double gsl_blas_dznrm2 (const gsl_vector_complex * X);
double gsl_blas_dzasum (const gsl_vector_complex * X);


CBLAS_INDEX_t gsl_blas_isamax (const gsl_vector_float * X);
CBLAS_INDEX_t gsl_blas_idamax (const gsl_vector * X);
CBLAS_INDEX_t gsl_blas_icamax (const gsl_vector_complex_float * X);
CBLAS_INDEX_t gsl_blas_izamax (const gsl_vector_complex * X);


int  gsl_blas_sswap (gsl_vector_float * X,
                     gsl_vector_float * Y);

int  gsl_blas_scopy (const gsl_vector_float * X,
                     gsl_vector_float * Y);

int  gsl_blas_saxpy (float alpha,
                     const gsl_vector_float * X,
                     gsl_vector_float * Y);

int  gsl_blas_dswap (gsl_vector * X,
                     gsl_vector * Y);

int  gsl_blas_dcopy (const gsl_vector * X,
                     gsl_vector * Y);

int  gsl_blas_daxpy (double alpha,
                     const gsl_vector * X,
                     gsl_vector * Y);

int  gsl_blas_cswap (gsl_vector_complex_float * X,
                     gsl_vector_complex_float * Y);

int  gsl_blas_ccopy (const gsl_vector_complex_float * X,
                     gsl_vector_complex_float * Y);

int  gsl_blas_caxpy (const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     gsl_vector_complex_float * Y);

int  gsl_blas_zswap (gsl_vector_complex * X,
                     gsl_vector_complex * Y);

int  gsl_blas_zcopy (const gsl_vector_complex * X,
                     gsl_vector_complex * Y);

int  gsl_blas_zaxpy (const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     gsl_vector_complex * Y);


int  gsl_blas_srotg (float a[], float b[], float c[], float s[]);

int  gsl_blas_srotmg (float d1[], float d2[], float b1[], float b2, float P[]);

int  gsl_blas_srot (gsl_vector_float * X,
                    gsl_vector_float * Y,
                    float c, float s);

int  gsl_blas_srotm (gsl_vector_float * X,
                     gsl_vector_float * Y,
                     const float P[]);

int  gsl_blas_drotg (double a[], double b[], double c[], double s[]);

int  gsl_blas_drotmg (double d1[], double d2[], double b1[],
                      double b2, double P[]);

int  gsl_blas_drot (gsl_vector * X,
                    gsl_vector * Y,
                    const double c, const double s);

int  gsl_blas_drotm (gsl_vector * X,
                     gsl_vector * Y,
                     const double P[]);


void gsl_blas_sscal  (float  alpha, gsl_vector_float * X);
void gsl_blas_dscal  (double alpha, gsl_vector * X);
void gsl_blas_cscal  (const gsl_complex_float * alpha, gsl_vector_complex_float * X);
void gsl_blas_zscal  (const gsl_complex * alpha, gsl_vector_complex * X);
void gsl_blas_csscal (float  alpha, gsl_vector_complex_float * X);
void gsl_blas_zdscal (double alpha, gsl_vector_complex * X);



/* ===========================================================================
 * Level 2
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
int  gsl_blas_sgemv (CBLAS_TRANSPOSE_t TransA,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y);

int  gsl_blas_sgbmv (CBLAS_TRANSPOSE_t TransA,
                     int KL, int KU,
		     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y);

int  gsl_blas_strmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X);

int  gsl_blas_stbmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X);

int  gsl_blas_stpmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const float Ap[],
                     gsl_vector_float * X);

int  gsl_blas_strsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X);

int  gsl_blas_stbsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X);

int  gsl_blas_stpsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const float Ap[],
                     gsl_vector_float * X);

int  gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y);

int  gsl_blas_dgbmv (CBLAS_TRANSPOSE_t TransA,
                     int KL, int KU,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y);

int  gsl_blas_dtrmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix * A,
                     gsl_vector * X);

int  gsl_blas_dtbmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix * A,
                     gsl_vector * X);

int  gsl_blas_dtpmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const double Ap[],
                     gsl_vector * X);

int  gsl_blas_dtrsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix * A,
                     gsl_vector * X);

int  gsl_blas_dtbsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix * A,
                     gsl_vector * X);

int  gsl_blas_dtpsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const double Ap[],
                     gsl_vector * X);

int  gsl_blas_cgemv (CBLAS_TRANSPOSE_t TransA,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y);

int  gsl_blas_cgbmv (CBLAS_TRANSPOSE_t TransA,
                     int KL, int KU,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y);

int  gsl_blas_ctrmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X);

int  gsl_blas_ctbmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X);

int  gsl_blas_ctpmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const void * Ap,
                     gsl_vector_complex_float * X);

int  gsl_blas_ctrsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X);

int  gsl_blas_ctbsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X);

int  gsl_blas_ctpsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const void * Ap,
                     gsl_vector_complex_float * X);

int  gsl_blas_zgemv (CBLAS_TRANSPOSE_t TransA,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y);

int  gsl_blas_zgbmv (CBLAS_TRANSPOSE_t TransA,
                     int KL, int KU,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y);

int  gsl_blas_ztrmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X);

int  gsl_blas_ztbmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X);

int  gsl_blas_ztpmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const void *Ap,
                     gsl_vector_complex * X);

int  gsl_blas_ztrsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex *X);

int  gsl_blas_ztbsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X);

int  gsl_blas_ztpsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const void *Ap,
                     gsl_vector_complex * X);

/*
 * Routines with S and D prefixes only
 */
int  gsl_blas_ssymv (CBLAS_UPLO_t Uplo,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y);

int  gsl_blas_ssbmv (CBLAS_UPLO_t Uplo,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y);

int  gsl_blas_sspmv (CBLAS_UPLO_t Uplo,
                     float alpha,
                     const float Ap[],
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y);

int  gsl_blas_sger (float alpha,
                    const gsl_vector_float * X,
                    const gsl_vector_float * Y,
                    gsl_matrix_float * A);

int  gsl_blas_ssyr (CBLAS_UPLO_t Uplo,
                    float alpha,
                    const gsl_vector_float * X,
                    gsl_matrix_float * A);

int  gsl_blas_sspr (CBLAS_UPLO_t Uplo,
                    float alpha,
                    const gsl_vector_float * X,
                    float Ap[]);

int  gsl_blas_ssyr2 (CBLAS_UPLO_t Uplo,
                     float alpha,
                     const gsl_vector_float * X,
                     const gsl_vector_float * Y,
                     gsl_matrix_float * A);

int  gsl_blas_sspr2 (CBLAS_UPLO_t Uplo,
                     float alpha,
                     const gsl_vector_float * X,
                     const gsl_vector_float * Y,
                     float * A);

int  gsl_blas_dsymv (CBLAS_UPLO_t Uplo,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y);

int  gsl_blas_dsbmv (CBLAS_UPLO_t Uplo,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y);

int  gsl_blas_dspmv (CBLAS_UPLO_t Uplo,
                     double alpha,
                     const double Ap[],
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y);

int  gsl_blas_dger (double alpha,
                    const gsl_vector * X,
                    const gsl_vector * Y,
                    gsl_matrix * A);

int  gsl_blas_dsyr (CBLAS_UPLO_t Uplo,
                    double alpha,
                    const gsl_vector * X,
                    gsl_matrix * A);

int  gsl_blas_dspr (CBLAS_UPLO_t Uplo,
                    double alpha,
                    const gsl_vector * X,
                    double Ap[]);

int  gsl_blas_dsyr2 (CBLAS_UPLO_t Uplo,
                     double alpha,
                     const gsl_vector * X,
                     const gsl_vector * Y,
                     gsl_matrix * A);

int  gsl_blas_dspr2 (CBLAS_UPLO_t Uplo,
                     double alpha,
                     const gsl_vector * X,
                     const gsl_vector * Y,
                     double * Ap);

/*
 * Routines with C and Z prefixes only
 */
int  gsl_blas_chemv (CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y);

int  gsl_blas_chbmv (CBLAS_UPLO_t Uplo,
                     int K,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y);

int  gsl_blas_chpmv (CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
		     const void * Ap,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y);

int  gsl_blas_cgeru (const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_matrix_complex_float * A);

int  gsl_blas_cgerc (const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_matrix_complex_float * A);

int  gsl_blas_cher (CBLAS_UPLO_t Uplo,
                    float alpha,
                    const gsl_vector_complex_float * X,
                    gsl_matrix_complex_float * A);

int  gsl_blas_chpr (CBLAS_UPLO_t Uplo,
                    float alpha,
                    const gsl_vector_complex_float * X,
                    void * Ap);

int  gsl_blas_cher2 (CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_matrix_complex_float * A);

int  gsl_blas_chpr2 (CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     void * Ap);

int  gsl_blas_zhemv (CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y);

int  gsl_blas_zhbmv (CBLAS_UPLO_t Uplo,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y);

int  gsl_blas_zhpmv (CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const void * Ap,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y);

int  gsl_blas_zgeru (const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A);

int  gsl_blas_zgerc (const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A);

int  gsl_blas_zher (CBLAS_UPLO_t Uplo,
                    double alpha,
                    const gsl_vector_complex * X,
                    gsl_matrix_complex * A);

int  gsl_blas_zhpr (CBLAS_UPLO_t Uplo,
                    double alpha,
                    const gsl_vector_complex * X,
                    void * Ap);

int  gsl_blas_zher2 (CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A);

int  gsl_blas_zhpr2 (CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     void *Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
int  gsl_blas_sgemm (CBLAS_TRANSPOSE_t TransA,
                     CBLAS_TRANSPOSE_t TransB,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_matrix_float * B,
                     float beta,
                     gsl_matrix_float * C);

int  gsl_blas_ssymm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_matrix_float * B,
                     float beta,
                     gsl_matrix_float * C);

int  gsl_blas_ssyrk (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     float beta,
                     gsl_matrix_float * C);

int  gsl_blas_ssyr2k (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans,
                      int K,
                      float alpha,
                      const gsl_matrix_float * A,
                      const gsl_matrix_float * B,
                      float beta,
                      gsl_matrix_float * C);

int  gsl_blas_strmm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     float alpha,
                     const gsl_matrix_float * A,
                     gsl_matrix_float * B);

int  gsl_blas_strsm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     float alpha,
                     const gsl_matrix_float * A,
                     gsl_matrix_float * B);

int  gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA,
                     CBLAS_TRANSPOSE_t TransB,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_matrix * B,
                     double beta,
                     gsl_matrix * C);

int  gsl_blas_dsymm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_matrix * B,
                     double beta,
                     gsl_matrix * C);

int  gsl_blas_dsyrk (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t Trans,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     double beta,
                     gsl_matrix * C);

int  gsl_blas_dsyr2k (CBLAS_UPLO_t Uplo,
                      CBLAS_TRANSPOSE_t Trans,
                      int K,
                      double alpha,
                      const  gsl_matrix * A,
                      const  gsl_matrix * B,
                      double beta,
                      gsl_matrix * C);

int  gsl_blas_dtrmm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     double alpha,
                     const gsl_matrix * A,
                     gsl_matrix * B);

int  gsl_blas_dtrsm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     double alpha,
                     const gsl_matrix * A,
                     gsl_matrix * B);

int  gsl_blas_cgemm (CBLAS_TRANSPOSE_t TransA,
                     CBLAS_TRANSPOSE_t TransB,
                     int K,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_matrix_complex_float * B,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C);

int  gsl_blas_csymm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_matrix_complex_float * B,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C);

int  gsl_blas_csyrk (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t Trans,
                     int K,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C);

int  gsl_blas_csyr2k (CBLAS_UPLO_t Uplo,
                      CBLAS_TRANSPOSE_t Trans,
                      int K,
                      const gsl_complex_float * alpha,
                      const gsl_matrix_complex_float * A,
                      const gsl_matrix_complex_float * B,
                      const gsl_complex_float * beta,
                      gsl_matrix_complex_float * C);

int  gsl_blas_ctrmm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     gsl_matrix_complex_float * B);

int  gsl_blas_ctrsm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     gsl_matrix_complex_float * B);

int  gsl_blas_zgemm (CBLAS_TRANSPOSE_t TransA,
                     CBLAS_TRANSPOSE_t TransB,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C);

int  gsl_blas_zsymm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C);

int  gsl_blas_zsyrk (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t Trans,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C);

int  gsl_blas_zsyr2k (CBLAS_UPLO_t Uplo,
                      CBLAS_TRANSPOSE_t Trans,
                      int K,
                      const gsl_complex * alpha,
                      const gsl_matrix_complex * A,
                      const gsl_matrix_complex * B,
                      const gsl_complex * beta,
                      gsl_matrix_complex *C);

int  gsl_blas_ztrmm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     gsl_matrix_complex * B);

int  gsl_blas_ztrsm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     gsl_matrix_complex * B);

/*
 * Routines with prefixes C and Z only
 */
int  gsl_blas_chemm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_matrix_complex_float * B,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C);

int  gsl_blas_cherk (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t Trans,
                     int K,
                     float alpha,
                     const gsl_matrix_complex_float * A,
                     float beta,
                     gsl_matrix_complex_float * C);

int  gsl_blas_cher2k (CBLAS_UPLO_t Uplo,
                      CBLAS_TRANSPOSE_t Trans,
                      int K,
                      const gsl_complex_float * alpha,
                      const gsl_matrix_complex_float * A,
                      const gsl_matrix_complex_float * B,
                      float beta,
                      gsl_matrix_complex_float * C);

int  gsl_blas_zhemm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C);

int  gsl_blas_zherk (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t Trans,
                     int K,
                     double alpha,
                     const gsl_matrix_complex * A,
                     double beta,
                     gsl_matrix_complex * C);

int  gsl_blas_zher2k (CBLAS_UPLO_t Uplo,
                      CBLAS_TRANSPOSE_t Trans,
                      int K,
                      const gsl_complex * alpha,
                      const gsl_matrix_complex * A,
                      const gsl_matrix_complex * B,
                      double beta,
                      gsl_matrix_complex * C);


#endif /* !GSL_BLAS_H_ */
