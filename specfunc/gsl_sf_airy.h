/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_AIRY_H_
#define GSL_SF_AIRY_H_


/* Airy function Ai(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_airy_Ai_impl(double x, double * result);
int     gsl_sf_airy_Ai_e(double x, double * result);
double  gsl_sf_airy_Ai(double x);


/* Airy function Bi(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int     gsl_sf_airy_Bi_impl(double x, double * result);
int     gsl_sf_airy_Bi_e(double x, double * result);
double  gsl_sf_airy_Bi(double x);


/* scaled Ai(x):
 *                     Ai(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai(x)   x > 0
 *
 * exceptions: none
 */
int     gsl_sf_airy_Ai_scaled_impl(double x, double * result);
int     gsl_sf_airy_Ai_scaled_e(double x, double * result);
double  gsl_sf_airy_Ai_scaled(double x);


/* scaled Bi(x):
 *                     Bi(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_scaled_impl(double x, double * result);
int gsl_sf_airy_Bi_scaled_e(double x, double * result);
double gsl_sf_airy_Bi_scaled(double x);


/* derivative Ai'(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_airy_Ai_deriv_impl(double x, double * result);
int     gsl_sf_airy_Ai_deriv_e(double x, double * result);
double  gsl_sf_airy_Ai_deriv(double x);


/* derivative Bi'(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int     gsl_sf_airy_Bi_deriv_impl(double x, double * result);
int     gsl_sf_airy_Bi_deriv_e(double x, double * result);
double  gsl_sf_airy_Bi_deriv(double x);


/* scaled derivative Ai'(x):
 *                     Ai'(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai'(x)   x > 0
 *
 * exceptions: none
 */
int     gsl_sf_airy_Ai_deriv_scaled_impl(double x, double * result);
int     gsl_sf_airy_Ai_deriv_scaled_e(double x, double * result);
double  gsl_sf_airy_Ai_deriv_scaled(double x);


/* scaled derivative:
 *                     Bi'(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi'(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_deriv_scaled_impl(double x, double * result);
int gsl_sf_airy_Bi_deriv_scaled_e(double x, double * result);
double gsl_sf_airy_Bi_deriv_scaled(double x);


#endif /* GSL_SF_AIRY_H_ */
