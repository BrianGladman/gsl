/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_LEGENDRE_H__
#define __GSL_SF_LEGENDRE_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


/* P_l(x)   l >= 0; |x| <= 1
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_legendre_Pl_impl(const int l, const double x, gsl_sf_result * result);
int     gsl_sf_legendre_Pl_e(const int l, const double x, gsl_sf_result * result);


/* P_l(x) for l=0,...,lmax; |x| <= 1
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_legendre_Pl_array_impl(const int lmax, const double x, double * result_array);
int gsl_sf_legendre_Pl_array_e(const int lmax, const double x, double * result_array);


/* P_l(x), l=1,2,3
 *
 * exceptions: none
 */
int gsl_sf_legendre_P1_impl(double x, gsl_sf_result * result);
int gsl_sf_legendre_P2_impl(double x, gsl_sf_result * result);
int gsl_sf_legendre_P3_impl(double x, gsl_sf_result * result);
int gsl_sf_legendre_P1_e(const double x, gsl_sf_result * result);
int gsl_sf_legendre_P2_e(const double x, gsl_sf_result * result);
int gsl_sf_legendre_P3_e(const double x, gsl_sf_result * result);


/* Q_0(x), x > -1, x != 1
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_legendre_Q0_impl(const double x, gsl_sf_result * result);
int gsl_sf_legendre_Q0_e(const double x, gsl_sf_result * result);


/* Q_1(x), x > -1, x != 1
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_legendre_Q1_impl(const double x, gsl_sf_result * result);
int gsl_sf_legendre_Q1_e(const double x, gsl_sf_result * result);


/* Q_l(x), x > -1, x != 1, l >= 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_legendre_Ql_impl(const int l, const double x, gsl_sf_result * result);
int gsl_sf_legendre_Ql_e(const int l, const double x, gsl_sf_result * result);


/* P_l^m(x)  m >= 0; l >= m; |x| <= 1.0
 *
 * Note that this function grows combinatorially with l.
 * Therefore we can easily generate an overflow for l larger
 * than about 150.
 *
 * There is no trouble for small m, but when m and l are both large,
 * then there will be trouble. Rather than allow overflows, these
 * functions refuse to calculate when they can sense that l and m are
 * too big.
 *
 * If you really want to calculate a spherical harmonic, then DO NOT
 * use this. Instead use legendre_sphPlm() below, which  uses a similar
 * recursion, but with the normalized functions.
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int     gsl_sf_legendre_Plm_impl(const int l, const int m, const double x, gsl_sf_result * result);
int     gsl_sf_legendre_Plm_e(const int l, const int m, const double x, gsl_sf_result * result);


/* P_l^m(x)  m >= 0; l >= m; |x| <= 1.0
 * l=|m|,...,lmax
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_legendre_Plm_array_impl(const int lmax, const int m, const double x, double * result_array);
int gsl_sf_legendre_Plm_array_e(const int lmax, const int m, const double x, double * result_array);


/* P_l^m(x), normalized properly for use in spherical harmonics
 * m >= 0; l >= m; |x| <= 1.0
 *
 * There is no overflow problem, as there is for the
 * standard normalization of P_l^m(x).
 *
 * Specifically, it returns:
 *
 *        sqrt((2l+1)/(4pi)) sqrt((l-m)!/(l+m)!) P_l^m(x)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_legendre_sphPlm_impl(const int l, int m, const double x, gsl_sf_result * result);
int     gsl_sf_legendre_sphPlm_e(const int l, const int m, const double x, gsl_sf_result * result);


/* P_l^m(x), normalized properly for use in spherical harmonics
 * m >= 0; l >= m; |x| <= 1.0
 * l=|m|,...,lmax
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_legendre_sphPlm_array_impl(const int lmax, int m, const double x, double * result_array);
int gsl_sf_legendre_sphPlm_array_e(const int lmax, const int m, const double x, double * result_array);


/* size of result_array[] needed for the array versions of Plm
 * (lmax - m + 1)
 */
int gsl_sf_legendre_array_size(const int lmax, const int m);


/* Irregular Spherical Conical Function
 * P^{1/2}_{-1/2 + I lambda}(x)
 *
 * x > -1.0
 * exceptions: GSL_EDOM
 */
int gsl_sf_conicalP_half_impl(const double lambda, const double x, gsl_sf_result * result);
int gsl_sf_conicalP_half_e(const double lambda, const double x, gsl_sf_result * result);


/* Regular Spherical Conical Function
 * P^{-1/2}_{-1/2 + I lambda}(x)
 *
 * x > -1.0
 * exceptions: GSL_EDOM
 */
int gsl_sf_conicalP_mhalf_impl(const double lambda, const double x, gsl_sf_result * result);
int gsl_sf_conicalP_mhalf_e(const double lambda, const double x, gsl_sf_result * result);


/* Conical Function
 * P^{0}_{-1/2 + I lambda}(x)
 *
 * x > -1.0
 * exceptions: GSL_EDOM
 */
int gsl_sf_conicalP_0_impl(const double lambda, const double x, gsl_sf_result * result);
int gsl_sf_conicalP_0_e(const double lambda, const double x, gsl_sf_result * result);


/* Conical Function
 * P^{1}_{-1/2 + I lambda}(x)
 *
 * x > -1.0
 * exceptions: GSL_EDOM
 */
int gsl_sf_conicalP_1_impl(const double lambda, const double x, gsl_sf_result * result);
int gsl_sf_conicalP_1_e(const double lambda, const double x, gsl_sf_result * result);


/* Regular Spherical Conical Function
 * P^{-1/2-l}_{-1/2 + I lambda}(x)
 *
 * x > -1.0, l >= -1
 * exceptions: GSL_EDOM
 */
int gsl_sf_conicalP_sph_reg_impl(const int l, const double lambda, const double x, gsl_sf_result * result);
int gsl_sf_conicalP_sph_reg_e(const int l, const double lambda, const double x, gsl_sf_result * result);


/* Regular Cylindrical Conical Function
 * P^{-m}_{-1/2 + I lambda}(x)
 *
 * x > -1.0, m >= -1
 * exceptions: GSL_EDOM
 */
int gsl_sf_conicalP_cyl_reg_impl(const int m, const double lambda, const double x, gsl_sf_result * result);
int gsl_sf_conicalP_cyl_reg_e(const int m, const double lambda, const double x, gsl_sf_result * result);


/* The following spherical functions are specializations
 * of Legendre functions which give the regular eigenfunctions
 * of the Laplacian on a 3-dimensional hyperbolic space.
 * Of particular interest is the flat limit, which is
 * Flat-Lim := {lambda->Inf, eta->0, lambda*eta fixed}.
 */
  
/* Zeroth radial eigenfunction of the Laplacian on the
 * 3-dimensional hyperbolic space.
 *
 * legendre_H3d_0(lambda,eta) := sin(lambda*eta)/(lambda*sinh(eta))
 * 
 * Normalization:
 * Flat-Lim legendre_H3d_0(lambda,eta) = j_0(lambda*eta)
 *
 * eta >= 0.0
 * exceptions: GSL_EDOM
 */
int gsl_sf_legendre_H3d_0_impl(const double lambda, const double eta, gsl_sf_result * result);
int gsl_sf_legendre_H3d_0_e(const double lambda, const double eta, gsl_sf_result * result);


/* First radial eigenfunction of the Laplacian on the
 * 3-dimensional hyperbolic space.
 *
 * legendre_H3d_1(lambda,eta) :=
 *    1/sqrt(lambda^2 + 1) sin(lam eta)/(lam sinh(eta))
 *    (coth(eta) - lambda cot(lambda*eta))
 * 
 * Normalization:
 * Flat-Lim legendre_H3d_1(lambda,eta) = j_1(lambda*eta)
 *
 * eta >= 0.0
 * exceptions: GSL_EDOM
 */
int gsl_sf_legendre_H3d_1_impl(const double lambda, const double eta, gsl_sf_result * result);
int gsl_sf_legendre_H3d_1_e(const double lambda, const double eta, gsl_sf_result * result);


/* l'th radial eigenfunction of the Laplacian on the
 * 3-dimensional hyperbolic space.
 *
 * Normalization:
 * Flat-Lim legendre_H3d_l(l,lambda,eta) = j_l(lambda*eta)
 *
 * eta >= 0.0, l >= 0
 * exceptions: GSL_EDOM
 */
int gsl_sf_legendre_H3d_impl(const int l, const double lambda, const double eta, gsl_sf_result * result);
int gsl_sf_legendre_H3d_e(const int l, const double lambda, const double eta, gsl_sf_result * result);


/* Array of H3d(ell),  0 <= ell <= lmax
 */
int gsl_sf_legendre_H3d_array_impl(const int lmax, const double lambda, const double eta, double * result_array);
int gsl_sf_legendre_H3d_array_e(const int lmax, const double lambda, const double eta, double * result_array);


#ifdef HAVE_INLINE
extern inline
int
gsl_sf_legendre_array_size(const int lmax, const int m)
{
  return lmax-m+1;
}
#endif /* HAVE_INLINE */


__END_DECLS

#endif /* __GSL_SF_LEGENDRE_H__ */
