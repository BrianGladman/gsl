/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_LEGENDRE_H_
#define GSL_LEGENDRE_H_

/* P_l(x)   l >= 0; x >= 0  */

int    gsl_sf_legendre_Pl_e(int l, double x, double * result);                /* GSL_EDOM */
int    gsl_sf_legendre_Pl_array_e(int lmax, double x, double * result_array); /* GSL_EDOM */ 

double gsl_sf_legendre_Pl(int l, double x);  /* domain */


/* P_l^m(x)  m >= 0; l >= m; x >= 0
 *
 * This function is a pain in the ass. The problem is that it grows
 * combinatorially with l. Therefore we can easily generate an
 * overflow for l larger than about 150.
 *
 * There is no trouble for small m, but when m and l are both large,
 * then there will be trouble. Rather than allow overflows, these
 * functions refuse to calculate when they can sense that l and m are
 * too big.
 *
 * If you really want to calculate a spherical harmonic, then DO NOT
 * use this. Instead use legendre_sphPlm() below, which  uses a similar
 * recursion, but with the normalized functions.
 */

int    gsl_sf_legendre_Plm_e(int l, int m, double x, double * result);                 /* GSL_EDOM, GSL_EOVRFLW */
int    gsl_sf_legendre_Plm_array_e(int lmax, int m, double x, double * result_array);  /* GSL_EDOM, GSL_EOVRFLW */

double gsl_sf_legendre_Plm(int l, int m, double x);  /* domain, overflow */


/* P_l^m(x), normalized properly for use in spherical harmonics
 * m >= 0; l >= m; x >= 0
 *
 * There is no overflow problem, as there is for the
 * standard normalization of P_l^m(x).
 *
 * Specifically, it returns:
 *
 *        sqrt((2l+1)/(4pi)) sqrt((l-m)!/(l+m)!) P_l^m(x)
 */
int    gsl_sf_legendre_sphPlm_e(int l, int m, double x, double * result);                /* GSL_EDOM */
int    gsl_sf_legendre_sphPlm_array_e(int lmax, int m, double x, double * result_array); /* GSL_EDOM */

double gsl_sf_legendre_sphPlm(int l, int m, double x);  /* domain */


/* size of result_array[] needed for the array versions (lmax - m + 1) */

int gsl_sf_legendre_array_size(int lmax, int m);


/* P_{-1/2 + I lambda}^{1/2}(x)  Irregular (Spherical) Conical Function */

int gsl_sf_conical_sph_irr_1_e(double lambda, double x, double * result);   /* GSL_EDOM */

double gsl_sf_conical_sph_irr_1(double lambda, double x);  /* domain */


/* P_{-1/2 + I lambda}^{-1/2}(x)  Regular (Spherical) Conical Function */

int gsl_sf_conical_sph_reg_0_e(double lambda, double x, double * result);   /* GSL_EDOM */

double gsl_sf_conical_sph_reg_0(double lambda, double x);  /* domain */


/* P_{-1/2 + I lambda}^{-1/2-l}(x)  Regular (Spherical) Conical Function */

int gsl_sf_conical_sph_reg_e(int l, double lambda, double x, double * result);        /* GSL_EDOM */
int gsl_sf_conical_sph_reg_array_e(int l, double lambda, double x, double * result);  /* GSL_EDOM */

double gsl_sf_conical_sph_reg(int l, double lambda, double x);  /* domain */


#endif /* !GSL_LEGENDRE_H_ */
