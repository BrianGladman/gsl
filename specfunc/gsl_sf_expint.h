/* Author: G. Jungman
 * RCS: $Id$
 */
#ifndef GSL_EXPINT_H_
#define GSL_EXPINT_H_


/* E_1(x)
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E1_impl(double x, double * result);
int     gsl_sf_expint_E1_e(double x, double * result);
double  gsl_sf_expint_E1(double x);


/* E_2(x)
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E2_impl(double x, double * result);
int     gsl_sf_expint_E2_e(double x, double * result);
double  gsl_sf_expint_E2(double x);


/* Ei(x)
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_Ei_impl(double x, double * result);
int     gsl_sf_expint_Ei_e(double x, double * result);
double  gsl_sf_expint_Ei(double x);


/* Shi(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_Shi_impl(double x, double * result);
int     gsl_sf_Shi_e(double x, double * result);
double  gsl_sf_Shi(double x);


/* Chi(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_Chi_impl(double x, double * result);
int     gsl_sf_Chi_e(double x, double * result);
double  gsl_sf_Chi(double x);


/* Si(x)
 *
 * exceptions: none
 */
int     gsl_sf_Si_impl(const double x, double * result);
int     gsl_sf_Si_e(double x, double * result);
double  gsl_sf_Si(double x);


/* Ci(x)
 *
 * exceptions: GSL_EDOM 
 */
int     gsl_sf_Ci_impl(const double x, double * result);
int     gsl_sf_Ci_e(double x, double * result);
double  gsl_sf_Ci(double x);



/* Ei_3(x) := Integral[ Exp[-t^3], {t,0,x}]
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_expint_3_impl(double x, double * result);
int     gsl_sf_expint_3_e(double x, double * result);
double  gsl_sf_expint_3(double x);


#endif /* !GSL_EXPINT_H_ */
