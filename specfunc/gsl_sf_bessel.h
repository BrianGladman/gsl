/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_BESSEL_H_
#define GSL_BESSEL_H_

/* Regular Bessel functions J_0(x), J_1(x), J_n(x) */

int gsl_sf_bessel_J0_e(double x, double * result);         /* none  */
int gsl_sf_bessel_J1_e(double x, double * result);         /* GSL_EUNDRFLW */
int gsl_sf_bessel_Jn_e(int n, double x, double * result);  /* GSL_EUNDRFLW */

double gsl_sf_bessel_J0(double x);          /* none  */
double gsl_sf_bessel_J1(double x);          /* underflow */
double gsl_sf_bessel_Jn(int n, double x);   /* underflow */


/* Regular modified Bessel functions I_0(x), I_1(x), I_n(x) */

int gsl_sf_bessel_I0_e(double x, double * result);        /* GSL_EOVRFLW */
int gsl_sf_bessel_I1_e(double x, double * result);        /* GSL_EOVRFLW, GSL_EUNDRFLW */
int gsl_sf_bessel_In_e(int n, double x, double * result); /* GSL_EOVRFLW, GSL_EUNDRFLW */

double gsl_sf_bessel_I0(double x);         /* overflow  */
double gsl_sf_bessel_I1(double x);         /* overflow, underflow */
double gsl_sf_bessel_In(int n, double x);  /* overflow, underflow */


/* Scaled regular modified Bessel funcions
 *  exp(-|x|) I_0(x)
 *  exp(-|x|) I_1(x)
 *  exp(-|x|) I_n(x)
 */

int gsl_sf_bessel_I0_scaled_e(double x, double * result);        /* none */
int gsl_sf_bessel_I1_scaled_e(double x, double * result);        /* GSL_EUNDRFLW */
int gsl_sf_bessel_In_scaled_e(int n, double x, double * result); /* GSL_EUNDRFLW */
int gsl_sf_bessel_In_scaled_array_e(int n, double x, double * result_array); /* GSL_EUNDRFLW */

double gsl_sf_bessel_I0_scaled(double x);        /* none       */
double gsl_sf_bessel_I1_scaled(double x);        /* underflow  */
double gsl_sf_bessel_In_scaled(int n, double x); /* underflow  */


/* Irregular Bessel functions Y_0(x), Y_1(x), Y_n(x)
 * x > 0.0
 */

int gsl_sf_bessel_Y0_e(double x, double * result);       /* GSL_EDOM, GSL_EUNDRFLW */
int gsl_sf_bessel_Y1_e(double x, double * result);       /* GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW */
int gsl_sf_bessel_Yn_e(int n,double x, double * result); /* GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW */

double gsl_sf_bessel_Y0(double x);         /* domain, underflow */
double gsl_sf_bessel_Y1(double x);         /* domain, overflow, underflow */
double gsl_sf_bessel_Yn(int n, double x);  /* domain, overflow, underflow */


/* Irregular modified Bessel functions K_0(x), K_1(x), K_n(x)
 *  x > 0.0
 */

int gsl_sf_bessel_K0_e(double x, double * result); /* GSL_EDOM, GSL_EUNDRFLW */
int gsl_sf_bessel_K1_e(double x, double * result); /* GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW */
int gsl_sf_bessel_Kn_e(int n, double x, double * result); 

double gsl_sf_bessel_K0(double x);          /* domain, underflow */
double gsl_sf_bessel_K1(double x);          /* domain, overflow, underflow */
double gsl_sf_bessel_Kn(int n, double x);


/* Scaled irregular modified Bessel functions
 *  exp(x) K_0(x)
 *  exp(x) K_1(x)
 *  exp(x) K_n(x)
 *  x > 0.0
 */

int gsl_sf_bessel_K0_scaled_e(double x, double * result);        /* GSL_EDOM */
int gsl_sf_bessel_K1_scaled_e(double x, double * result);        /* GSL_EDOM, GSL_EUNDRFLW */
int gsl_sf_bessel_Kn_scaled_e(int n, double x, double * result); /* GSL_EDOM, GSL_EUNDRFLW */

double gsl_sf_bessel_K0_scaled(double x);        /* domain */
double gsl_sf_bessel_K1_scaled(double x);        /* domain, underflow */
double gsl_sf_bessel_Kn_scaled(int n, double x); /* domain, underflow */


/* Regular spherical Bessel functions j_0(x), j_1(x), j_2(x) */

int gsl_sf_bessel_j0_e(double x, double * result);   /* none */
int gsl_sf_bessel_j1_e(double x, double * result);   /* GSL_EUNDRFLW */
int gsl_sf_bessel_j2_e(double x, double * result);   /* GSL_EUNDRFLW */

double gsl_sf_bessel_j0(double x);   /* none      */
double gsl_sf_bessel_j1(double x);   /* underflow */
double gsl_sf_bessel_j2(double x);   /* underflow */


/* Regular spherical Bessel functions j_l(x) */

int gsl_sf_bessel_jl_e(int l, double x, double * result);
int gsl_sf_bessel_jl_array_e(int lmax, double x, double * result_array);

double gsl_sf_bessel_jl(int l, double x);


/* Irregular spherical Bessel functions y_0(x), y_1(x), y_2(x) */

int gsl_sf_bessel_y0_e(double x, double * result);   /* none */
int gsl_sf_bessel_y1_e(double x, double * result);   /* GSL_EUNDRFLW */
int gsl_sf_bessel_y2_e(double x, double * result);   /* GSL_EUNDRFLW */

double gsl_sf_bessel_y0(double x);   /* none      */
double gsl_sf_bessel_y1(double x);   /* underflow */
double gsl_sf_bessel_y2(double x);   /* underflow */


/* Irregular spherical Bessel functions y_l(x) */

int gsl_sf_bessel_yl_e(int l, double x, double * result);
int gsl_sf_bessel_yl_array_e(int lmax, double x, double * result_array);

double gsl_sf_bessel_yl(int l, double x);



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




int gsl_sf_bessel_I0_scaled_impl(double x, double * result);
int gsl_sf_bessel_I0_impl(double x, double * result);

int gsl_sf_bessel_I1_scaled_impl(double x, double * result);
int gsl_sf_bessel_I1_impl(double x, double * result);

int gsl_sf_bessel_In_scaled_impl(int n, double x, double * result);
int gsl_sf_bessel_In_impl(int n, double x, double * result);

int gsl_sf_bessel_In_scaled_array_impl(int nmax, double x, double * result_array);
int gsl_sf_bessel_In_array_impl(int nmax, double x, double * result_array);


int gsl_sf_bessel_J0_impl(double x, double * result);
int gsl_sf_bessel_J1_impl(double x, double * result);
int gsl_sf_bessel_Jn_impl(int n, double x, double * result);

int gsl_sf_bessel_K0_scaled_impl(double x, double * result);
int gsl_sf_bessel_K0_impl(double x, double * result);

int gsl_sf_bessel_K1_scaled_impl(double x, double * result);
int gsl_sf_bessel_K1_impl(double x, double * result);

int gsl_sf_bessel_Kn_scaled_impl(int n, double x, double * result);
int gsl_sf_bessel_Kn_impl(int n, double x, double * result);

int gsl_sf_bessel_Kn_scaled_array_impl(int nmax, double x, double * result_array);
int gsl_sf_bessel_Kn_array_impl(int nmax, double x, double * result_array);

int gsl_sf_bessel_Y0_impl(double x, double * result);
int gsl_sf_bessel_Y1_impl(double x, double * result);
int gsl_sf_bessel_Yn_impl(int n, double x, double * result);


int gsl_sf_bessel_j0_impl(double x, double * result);
int gsl_sf_bessel_j1_impl(double x, double * result);
int gsl_sf_bessel_j2_impl(double x, double * result);
int gsl_sf_bessel_jl_impl(int l, double x, double * result);
int gsl_sf_bessel_jl_array_impl(int lmax, double x, double * result_array);
int gsl_sf_bessel_j_steed_array_impl(int lmax, double x, double * jl_x);

int gsl_sf_bessel_y0_impl(double x, double * result);
int gsl_sf_bessel_y1_impl(double x, double * result);
int gsl_sf_bessel_y2_impl(double x, double * result);

int gsl_sf_bessel_yl_array_impl(int lmax, double x, double * result_array);


#endif /* !GSL_BESSEL_H_ */
