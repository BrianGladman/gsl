/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_BESSEL_H_
#define GSL_BESSEL_H_


/* Regular Bessel function J_0(x)
 *
 * exceptions: none
 */
int     gsl_sf_bessel_J0_impl(double x, double * result);
int     gsl_sf_bessel_J0_e(double x, double * result);
double  gsl_sf_bessel_J0(double x);


/* Regular Bessel function J_1(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_J1_impl(double x, double * result);
int     gsl_sf_bessel_J1_e(double x, double * result);
double  gsl_sf_bessel_J1(double x);


/* Regular Bessel function J_n(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_Jn_impl(int n, double x, double * result);
int     gsl_sf_bessel_Jn_e(int n, double x, double * result);
double  gsl_sf_bessel_Jn(int n, double x);


/* Irregular Bessel function Y_0(x)
 *
 * x > 0.0
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_Y0_impl(double x, double * result);
int     gsl_sf_bessel_Y0_e(double x, double * result);
double  gsl_sf_bessel_Y0(double x);


/* Irregular Bessel function Y_1(x)
 *
 * x > 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_Y1_impl(double x, double * result);
int     gsl_sf_bessel_Y1_e(double x, double * result);
double  gsl_sf_bessel_Y1(double x);


/* Irregular Bessel function Y_n(x)
 *
 * x > 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_Yn_impl(int n,double x, double * result);
int     gsl_sf_bessel_Yn_e(int n,double x, double * result);
double  gsl_sf_bessel_Yn(int n, double x);


/* Regular modified Bessel function I_0(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int     gsl_sf_bessel_I0_impl(double x, double * result);
int     gsl_sf_bessel_I0_e(double x, double * result);
double  gsl_sf_bessel_I0(double x);


/* Regular modified Bessel function I_1(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_I1_impl(double x, double * result);
int     gsl_sf_bessel_I1_e(double x, double * result);
double  gsl_sf_bessel_I1(double x);


/* Regular modified Bessel function I_n(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_In_impl(int n, double x, double * result);
int     gsl_sf_bessel_In_e(int n, double x, double * result);
double  gsl_sf_bessel_In(int n, double x);


/* Scaled regular modified Bessel function
 *  exp(-|x|) I_0(x)
 *
 * exceptions: none
 */
int     gsl_sf_bessel_I0_scaled_impl(double x, double * result);
int     gsl_sf_bessel_I0_scaled_e(double x, double * result);
double  gsl_sf_bessel_I0_scaled(double x);


/* Scaled regular modified Bessel function
 *  exp(-|x|) I_1(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_I1_scaled_impl(double x, double * result);
int     gsl_sf_bessel_I1_scaled_e(double x, double * result);
double  gsl_sf_bessel_I1_scaled(double x);


/* Scaled regular modified Bessel function
 *  exp(-|x|) I_n(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_In_scaled_impl(int n, double x, double * result);
int     gsl_sf_bessel_In_scaled_e(int n, double x, double * result);
double  gsl_sf_bessel_In_scaled(int n, double x);


/* Scaled regular modified Bessel function
 *  exp(-|x|) I_k(x)  for k=0,1,...,n
 *
 * exceptions: GSL_EUNDRFLW
 */
int  gsl_sf_bessel_In_scaled_array_impl(int n, double x, double * result_array);
int  gsl_sf_bessel_In_scaled_array_e(int n, double x, double * result_array);


/* Irregular modified Bessel function K_0(x)
 *
 * x > 0.0
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_K0_impl(double x, double * result);
int     gsl_sf_bessel_K0_e(double x, double * result);
double  gsl_sf_bessel_K0(double x);


/* Irregular modified Bessel function K_1(x)
 *
 * x > 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_K1_impl(double x, double * result);
int     gsl_sf_bessel_K1_e(double x, double * result);
double  gsl_sf_bessel_K1(double x);


/* Irregular modified Bessel function K_n(x)
 *
 * x > 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_Kn_impl(int n, double x, double * result);
int     gsl_sf_bessel_Kn_e(int n, double x, double * result);
double  gsl_sf_bessel_Kn(int n, double x);


/* Scaled irregular modified Bessel function
 *  exp(x) K_0(x)
 *
 * x > 0.0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_bessel_K0_scaled_impl(double x, double * result);
int     gsl_sf_bessel_K0_scaled_e(double x, double * result);
double  gsl_sf_bessel_K0_scaled(double x);


/* Scaled irregular modified Bessel function
 *  exp(x) K_1(x)
 *
 * x > 0.0
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_K1_scaled_impl(double x, double * result); 
int     gsl_sf_bessel_K1_scaled_e(double x, double * result);
double  gsl_sf_bessel_K1_scaled(double x);


/* Scaled irregular modified Bessel function
 *  exp(x) K_n(x)
 *
 * x > 0.0
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_Kn_scaled_impl(int n, double x, double * result);
int     gsl_sf_bessel_Kn_scaled_e(int n, double x, double * result);
double  gsl_sf_bessel_Kn_scaled(int n, double x);


/* Regular spherical Bessel function j_0(x)
 *
 * exceptions: none
 */
int     gsl_sf_bessel_j0_impl(double x, double * result);
int     gsl_sf_bessel_j0_e(double x, double * result);
double  gsl_sf_bessel_j0(double x);


/* Regular spherical Bessel function j_1(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_j1_impl(double x, double * result);
int     gsl_sf_bessel_j1_e(double x, double * result);
double  gsl_sf_bessel_j1(double x);


/* Regular spherical Bessel function j_2(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_j2_impl(double x, double * result);
int     gsl_sf_bessel_j2_e(double x, double * result);
double  gsl_sf_bessel_j2(double x);


/* Regular spherical Bessel function j_l(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_jl_impl(int l, double x, double * result);
int     gsl_sf_bessel_jl_e(int l, double x, double * result);
double  gsl_sf_bessel_jl(int l, double x);


/* Regular spherical Bessel function j_l(x) for l=0,1,...,lmax
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_bessel_jl_array_impl(int lmax, double x, double * result_array);
int gsl_sf_bessel_jl_array_e(int lmax, double x, double * result_array);
int gsl_sf_bessel_j_steed_array_impl(int lmax, double x, double * jl_x);


/* Irregular spherical Bessel function y_0(x)
 *
 * exceptions: none
 */
int     gsl_sf_bessel_y0_impl(double x, double * result);
int     gsl_sf_bessel_y0_e(double x, double * result);
double  gsl_sf_bessel_y0(double x);


/* Irregular spherical Bessel function y_1(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_y1_impl(double x, double * result);
int     gsl_sf_bessel_y1_e(double x, double * result);
double  gsl_sf_bessel_y1(double x);


/* Irregular spherical Bessel function y_2(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_y2_impl(double x, double * result);
int     gsl_sf_bessel_y2_e(double x, double * result);
double  gsl_sf_bessel_y2(double x);


/* Irregular spherical Bessel function y_l(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_yl_impl(int l, double x, double * result);
int     gsl_sf_bessel_yl_e(int l, double x, double * result);
double  gsl_sf_bessel_yl(int l, double x);


/* Irregular spherical Bessel function y_l(x) for l=0,1,...,lmax
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_bessel_yl_array_impl(int lmax, double x, double * result_array);
int gsl_sf_bessel_yl_array_e(int lmax, double x, double * result_array);


/* Regular scaled modified spherical Bessel function
 *
 * Exp[-|x|] i_0(x)
 *
 * exceptions: none
 */
int     gsl_sf_bessel_i0_scaled_impl(double x, double * result);
int     gsl_sf_bessel_i0_scaled_e(double x, double * result);
double  gsl_sf_bessel_i0_scaled(double x);


/* Regular scaled modified spherical Bessel function
 *
 * Exp[-|x|] i_1(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_i1_scaled_impl(double x, double * result);
int     gsl_sf_bessel_i1_scaled_e(double x, double * result);
double  gsl_sf_bessel_i1_scaled(double x);


/* Regular scaled modified spherical Bessel function
 *
 * Exp[-|x|] i_2(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_i2_scaled_impl(double x, double * result);
int     gsl_sf_bessel_i2_scaled_e(double x, double * result);
double  gsl_sf_bessel_i2_scaled(double x);


/* Regular scaled modified spherical Bessel functions
 *
 * Exp[-|x|] i_l(x)
 *
 * i_l(x) = Sqrt[Pi/(2x)] BesselI[l+1/2,x]
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_il_scaled_impl(int l, double x, double * result);
int     gsl_sf_bessel_il_scaled_e(int l, double x, double * result);
double  gsl_sf_bessel_il_scaled(int l, double x);


/* Regular scaled modified spherical Bessel functions
 *
 * i_l(x)
 * for l=0,1,...,lmax
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_bessel_il_scaled_array_impl(int lmax, double x, double * result_array);
int gsl_sf_bessel_il_scaled_array_e(int lmax, double x, double * result_array);


/* Irregular modified spherical Bessel function k_0(x)
 *
 * exceptions: none
 */
int     gsl_sf_bessel_k0_scaled_impl(double x, double * result);
int     gsl_sf_bessel_k0_scaled_e(double x, double * result);
double  gsl_sf_bessel_k0_scaled(double x);


/* Irregular modified spherical Bessel function k_1(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_k1_scaled_impl(double x, double * result);
int     gsl_sf_bessel_k1_scaled_e(double x, double * result);
double  gsl_sf_bessel_k1_scaled(double x);


/* Irregular modified spherical Bessel function k_2(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_k2_scaled_impl(double x, double * result);
int     gsl_sf_bessel_k2_scaled_e(double x, double * result);
double  gsl_sf_bessel_k2_scaled(double x);


/* Irregular modified spherical Bessel function
 *
 * k_l(x) = Sqrt[Pi/(2x)] BesselK[l+1/2,x]
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_bessel_kl_scaled_impl(int l, double x, double * result);
int     gsl_sf_bessel_kl_scaled_e(int l, double x, double * result);
double  gsl_sf_bessel_kl_scaled(int l, double x);


/* Irregular modified spherical Bessel function
 *
 * k_l(x)
 * for l=0,1,...,lmax
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_bessel_kl_scaled_array_impl(int lmax, double x, double * result_array);
int gsl_sf_bessel_kl_scaled_array_e(int lmax, double x, double * result_array);


/* Regular cylindrical Bessel function J_nu(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_Jnu_impl(double nu, double x, double * Jnu);
int     gsl_sf_bessel_Jnu_e(double nu, double x, double * Jnu);
double  gsl_sf_bessel_Jnu(double nu, double x);


/* Irregular cylindrical Bessel function Y_nu(x)
 *
 * exceptions:  
 */
int     gsl_sf_bessel_Ynu_impl(double nu, double x, double * result);
int     gsl_sf_bessel_Ynu_e(double nu, double x, double * result);
double  gsl_sf_bessel_Ynu(double nu, double x);


/* Scaled modified cylindrical Bessel functions
 *
 * Exp[-|x|] BesselI[nu, x]
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_bessel_Inu_scaled_impl(double nu, double x, double * result);
int     gsl_sf_bessel_Inu_scaled_e(double nu, double x, double * result);
double  gsl_sf_bessel_Inu_scaled(double nu, double x);


/* Modified cylindrical Bessel functions
 *
 * BesselI[nu, x]
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int     gsl_sf_bessel_Inu_impl(double nu, double x, double * result);
int     gsl_sf_bessel_Inu_e(double nu, double x, double * result);
double  gsl_sf_bessel_Inu(double nu, double x);



/* Scaled modified cylindrical Bessel functions
 *
 * Exp[+|x|] BesselK[nu, x]
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_bessel_Knu_scaled_impl(double nu, double x, double * result);
int     gsl_sf_bessel_Knu_scaled_e(double nu, double x, double * result);
double  gsl_sf_bessel_Knu_scaled(double nu, double x);


/* Modified cylindrical Bessel functions
 *
 * BesselK[nu, x]
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_bessel_Knu_impl(double nu, double x, double * result);
int     gsl_sf_bessel_Knu_e(double nu, double x, double * result);
double  gsl_sf_bessel_Knu(double nu, double x);


/* Regular cylindrical Bessel functions J_nu(x) calculated
   with the Meissel uniform approximation. Assumes x >= 0.
   This approximation is accurate to near 10^{-3} at the boundaries
   between the asymptotic regions; well away from the boundaries
   the accuracy is better than 10^{-5}.
   Asymptotic approximations 8.11.5, 8.12.5, and 8.42.7 from
   G.N.Watson, A Treatise on the Theory of Bessel Functions,
   2nd Edition (Cambridge University Press, 1944).
   Higher terms in expansion for x near l given by
   Airey in Phil. Mag. 31, 520 (1916).
 */
double besselJ_meissel(double nu, double x);


/* Evaluate  d/dx(J_nu(x)) given the value J_nu(x).
   Uses the recursion relation and one extra function evaluation.
   Assumes x > 0.
   This is not very good pointwise (although the other relation,
   written in terms of nu+1 and nu-1, is only slightly better).
   However, it is quite good in the mean, over intervals greater
   than one period.
 */
double besselJprime_meissel(double nu, double x, double J_nu_x);


/* Evaluate Meissel expansion for spherical Bessel functions j_l(x).
 */
double sphbesselj_meissel(double l, double x);


/* Evaluate d/dx(j_l(x)) using the value j_l(x). Use recursion and 
 * one extra function evaluation.
 */
double sphbesseljprime_meissel(double l, double x, double jl_x);


/* Leading large argument behaviour of the regular
   spherical Bessel function j_l(x), using the 
   Meissel approximation.
   Calculates the derivative as well if dflag!=0.
 */
void asymp_sphbesselj_meissel(double l, double x,
				     double *jl, double *jlp,
				     int dflag);


#endif /* !GSL_BESSEL_H_ */
