/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_AIRY_H_
#define GSL_SF_AIRY_H_


/* Airy function Ai(x) */

int gsl_sf_airy_Ai_e(double x, double * result);   /* GSL_EUNDRFLW */

double gsl_sf_airy_Ai(double x);   /* underflow */


/* Airy function Bi(x) */

int gsl_sf_airy_Bi_e(double x, double * result);   /* GSL_EOVRFLW */

double gsl_sf_airy_Bi(double x);  /* overflow */


/* scaled Ai(x):
 *                    Ai(x)   x < 0
 *  exp(+2/3 x^{3/2}) Ai(x)   x > 0
 */

int gsl_sf_airy_Ai_scaled_e(double x, double * result);  /* none */

double gsl_sf_airy_Ai_scaled(double x);  /* none */


/* scaled Bi(x):
 *                    Bi(x)   x < 0
 *  exp(-2/3 x^{3/2}) Bi(x)   x > 0
 */

int gsl_sf_airy_Bi_scaled_e(double x, double * result);  /* none */

double gsl_sf_airy_Bi_scaled(double x);  /* none */


/* derivative Ai'(x) */

int gsl_sf_airy_Ai_deriv_e(double x, double * result);  /* GSL_EUNDRFLW */

double gsl_sf_airy_Ai_deriv(double x);  /* underflow */


/* derivative Bi'(x) */

int gsl_sf_airy_Bi_deriv_e(double x, double * result);  /* GSL_EOVRFLW */

double gsl_sf_airy_Bi_deriv(double x);  /* overflow */


/* scaled derivative Ai'(x):
 *                    Ai'(x)   x < 0
 *  exp(+2/3 x^{3/2}) Ai'(x)   x > 0
 */

int gsl_sf_airy_Ai_deriv_scaled_e(double x, double * result);  /* none */

double gsl_sf_airy_Ai_deriv_scaled(double x);  /* none */


/* scaled derivative:
 *                    Bi'(x)   x < 0
 *  exp(-2/3 x^{3/2}) Bi'(x)   x > 0
 */

int gsl_sf_airy_Bi_deriv_scaled_e(double x, double * result);  /* none */

double gsl_sf_airy_Bi_deriv_scaled(double x);  /* none */



int gsl_sf_airy_Ai_impl(double x, double * result);
int gsl_sf_airy_Bi_impl(double x, double * result);
int gsl_sf_airy_Ai_scaled_impl(double x, double * result);
int gsl_sf_airy_Bi_scaled_impl(double x, double * result);

int gsl_sf_airy_Ai_deriv_scaled_impl(double x, double * result);
int gsl_sf_airy_Bi_deriv_scaled_impl(double x, double * result);
int gsl_sf_airy_Ai_deriv_impl(double x, double * result);
int gsl_sf_airy_Bi_deriv_impl(double x, double * result);


#endif /* GSL_SF_AIRY_H_ */
