/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_LEGENDRE_H_
#define GSL_LEGENDRE_H_


/* Calculate P_l^m(x).
 * Requires m positive, which is no loss of generality.
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
int    gsl_sf_legendre_Plm_e(int l, int m, double x, double * result);
double gsl_sf_legendre_Plm(int l, int m, double x);
int    gsl_sf_legendre_Plm_array_e(int lmax, int m, double x, double * result_array);


/* Calculate associated Legendre function P_l^m(x), but
 * normalized properly for use in spherical harmonics.
 * There is no overflow problem, as there is for the
 * standard normalization of P_l^m(x).
 *
 * Specifically, it returns:
 *
 *        sqrt((2l+1)/(4pi)) sqrt((l-m)!/(l+m)!) P_l^m(x)
 */
int    gsl_sf_legendre_sphPlm_e(int l, int m, double x, double * result);
double gsl_sf_legendre_sphPlm(int l, int m, double x);
int    gsl_sf_legendre_sphPlm_array_e(int lmax, int m, double x, double * result_array);


/* size of result_array[] needed for the array versions (lmax - m + 1) */
int gsl_sf_legendre_array_size(int lmax, int m);


#endif /* !GSL_LEGENDRE_H_ */
