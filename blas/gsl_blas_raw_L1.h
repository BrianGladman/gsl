/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
/* Prototypes for level 1 BLAS functions
 * Based on draft BLAST C interface specification  [Jul 7 1998]
 */
#ifndef GSL_BLAS_RAW_L1_H_
#define GSL_BLAS_RAW_L1_H_

#include <gsl_blas_types.h>


float  gsl_blas_raw_sdsdot (size_t N,
                            float alpha,
                            const float X[], int incX,
                            const float Y[], int incY);

double gsl_blas_raw_dsdot (size_t N,
                           const float X[], int incX,
                           const float Y[], int incY);

float  gsl_blas_raw_sdot (size_t N,
                          const float X[], int incX,
                          const float Y[], int incY);
                          
double gsl_blas_raw_ddot (size_t N,
                          const double X[], int incX,
                          const double Y[], int incY);

void gsl_blas_raw_cdotu (size_t N,
                         const gsl_complex_packed_array_float X, int incX,
                         const gsl_complex_packed_array_float Y, int incY,
		         gsl_complex_packed_float dotu);

void gsl_blas_raw_cdotc (size_t N,
                         const gsl_complex_packed_array_float X, int incX,
                         const gsl_complex_packed_array_float Y, int incY,
                         gsl_complex_packed_float dotc);

void gsl_blas_raw_zdotu (size_t N,
                         const gsl_complex_packed_array X, int incX,
                         const gsl_complex_packed_array Y, int incY,
                         gsl_complex_packed dotu);

void gsl_blas_raw_zdotc (size_t N,
                         const gsl_complex_packed_array X, int incX,
                         const gsl_complex_packed_array Y, int incY,
                         gsl_complex_packed dotc);

float  gsl_blas_raw_snrm2  (size_t N, const float  X[], int incX);
double gsl_blas_raw_dnrm2  (size_t N, const double X[], int incX);
float  gsl_blas_raw_scnrm2 (size_t N, const gsl_complex_packed_array_float X, int incX);
double gsl_blas_raw_dznrm2 (size_t N, const gsl_complex_packed_array X, int incX);
float  gsl_blas_raw_sasum  (size_t N, const float  X[], int incX);
double gsl_blas_raw_dasum  (size_t N, const double X[], int incX);
float  gsl_blas_raw_scasum (size_t N, const gsl_complex_packed_array_float X, int incX);
double gsl_blas_raw_dzasum (size_t N, const gsl_complex_packed_array X, int incX);

CBLAS_INDEX gsl_blas_raw_isamax (size_t N, const float  X[], int incX);
CBLAS_INDEX gsl_blas_raw_idamax (size_t N, const double X[], int incX);
CBLAS_INDEX gsl_blas_raw_icamax (size_t N, const gsl_complex_packed_array_float X, int incX);
CBLAS_INDEX gsl_blas_raw_izamax (size_t N, const gsl_complex_packed_array X, int incX);


void gsl_blas_raw_sswap (size_t N,
                         float X[], int incX,
                         float Y[], int incY);

void gsl_blas_raw_dswap (size_t N,
                         double X[], int incX,
                         double Y[], int incY);

void gsl_blas_raw_cswap (size_t N,
                         gsl_complex_packed_array_float X, int incX,
                         gsl_complex_packed_array_float Y, int incY);

void gsl_blas_raw_zswap (size_t N,
                         gsl_complex_packed_array X, int incX,
                         gsl_complex_packed_array Y, int incY);

void gsl_blas_raw_scopy (size_t N,
                         const float X[], int incX,
                         float Y[], int incY);

void gsl_blas_raw_dcopy (size_t N,
                         const double X[], int incX,
                         double Y[], int incY);

void gsl_blas_raw_ccopy (size_t N,
                         const gsl_complex_packed_array_float X, int incX,
                         gsl_complex_packed_array_float Y, int incY);

void gsl_blas_raw_zcopy (size_t N,
                         const gsl_complex_packed_array X, int incX,
                         gsl_complex_packed_array Y, int incY);

void gsl_blas_raw_saxpy (size_t N,
                         float alpha,
                         const float X[], int incX,
                         float Y[], int incY);


void gsl_blas_raw_daxpy (size_t N,
                         double alpha,
                         const double X[], int incX, 
                         double Y[], int incY);


void gsl_blas_raw_caxpy (size_t N,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float X, int incX,
                         gsl_complex_packed_array_float Y, int incY);

void gsl_blas_raw_zaxpy (size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array X, int incX,
                         gsl_complex_packed_array Y, int incY);

void gsl_blas_raw_srotg (float a[], float b[], float c[], float s[]);
void gsl_blas_raw_drotg (double a[], double b[], double c[], double s[]);

void gsl_blas_raw_srotmg (float d1[], float d2[], float b1[], float b2, float P[]);
void gsl_blas_raw_drotmg (double d1[], double d2[], double b1[],
                          double b2, double P[]);


void gsl_blas_raw_srot (size_t N,
                        float X[], int incX,
                        float Y[], int incY,
                        float c, float s);
void gsl_blas_raw_drot (size_t N,
                        double X[], int incX,
                        double Y[], int incY,
                        const double c, const double s);

void gsl_blas_raw_srotm (size_t N,
                         float X[], int incX,
                         float Y[], int incY,
                         const float P[]);
void gsl_blas_raw_drotm (size_t N,
                         double X[], int incX,
                         double Y[], int incY,
			 const double P[]);

void gsl_blas_raw_sscal  (size_t N, float  alpha, float  X[], int incX);
void gsl_blas_raw_dscal  (size_t N, double alpha, double X[], int incX);
void gsl_blas_raw_cscal  (size_t N, const gsl_complex_packed_float alpha, gsl_complex_packed_array_float X, int incX);
void gsl_blas_raw_zscal  (size_t N, const gsl_complex_packed alpha, gsl_complex_packed_array X, int incX);
void gsl_blas_raw_csscal (size_t N, float  alpha, gsl_complex_packed_array_float X, int incX);
void gsl_blas_raw_zdscal (size_t N, double alpha, gsl_complex_packed_array X, int incX);

#endif  /* !GSL_BLAS_RAW_L1_H_ */
