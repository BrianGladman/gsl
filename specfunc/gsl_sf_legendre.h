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
 * then there will be trouble. Rather than allow overflows, this
 * function refuses to calculate when it can sense that l and m are
 * too big. A warning is generated.
 *
 * If you really want to calculate a spherical harmonic, then DO NOT
 * use this. Instead use Ylm_legendre() below, which  uses a similar
 * recursion, but with the normalized functions.
 */
double gsl_sf_legendre(int l, int m, double x);


/* Calculate associated Legendre function P_l^m(x), but
 * normalized properly for use in spherical harmonics. This also
 * has the advantage that this function is not growing wildly with l,
 * so there is no overflow problem.
 *
 * Specifically, it returns:
 *
 *        sqrt((2l+1)/(4pi)) sqrt((l-m)!/(l+m)!) P_l^m(x)
 */
double gsl_sf_Ylm_legendre(int l, int m, double x);


#endif /* !GSL_LEGENDRE_H_ */
