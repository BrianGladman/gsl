/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_AIRY_H_
#define GSL_SF_AIRY_H_

#include "gsl_precision.h"


/* Airy function Ai(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_airy_Ai_impl(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
int     gsl_sf_airy_Ai_e(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
double  gsl_sf_airy_Ai(double x, gsl_prec_t goal, unsigned int err_bits);


/* Airy function Bi(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int     gsl_sf_airy_Bi_impl(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
int     gsl_sf_airy_Bi_e(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
double  gsl_sf_airy_Bi(double x, gsl_prec_t goal, unsigned int err_bits);


/* scaled Ai(x):
 *                     Ai(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai(x)   x > 0
 *
 * exceptions: none
 */
int     gsl_sf_airy_Ai_scaled_impl(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
int     gsl_sf_airy_Ai_scaled_e(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
double  gsl_sf_airy_Ai_scaled(double x, gsl_prec_t goal, unsigned int err_bits);


/* scaled Bi(x):
 *                     Bi(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_scaled_impl(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
int gsl_sf_airy_Bi_scaled_e(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
double gsl_sf_airy_Bi_scaled(double x, gsl_prec_t goal, unsigned int err_bits);


/* derivative Ai'(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_airy_Ai_deriv_impl(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
int     gsl_sf_airy_Ai_deriv_e(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
double  gsl_sf_airy_Ai_deriv(double x, gsl_prec_t goal, unsigned int err_bits);


/* derivative Bi'(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int     gsl_sf_airy_Bi_deriv_impl(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
int     gsl_sf_airy_Bi_deriv_e(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
double  gsl_sf_airy_Bi_deriv(double x, gsl_prec_t goal, unsigned int err_bits);


/* scaled derivative Ai'(x):
 *                     Ai'(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai'(x)   x > 0
 *
 * exceptions: none
 */
int     gsl_sf_airy_Ai_deriv_scaled_impl(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
int     gsl_sf_airy_Ai_deriv_scaled_e(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
double  gsl_sf_airy_Ai_deriv_scaled(double x, gsl_prec_t goal, unsigned int err_bits);


/* scaled derivative:
 *                     Bi'(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi'(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_deriv_scaled_impl(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
int gsl_sf_airy_Bi_deriv_scaled_e(double x, double * result, gsl_prec_t goal, unsigned int err_bits);
double gsl_sf_airy_Bi_deriv_scaled(double x, gsl_prec_t goal, unsigned int err_bits);


#endif /* GSL_SF_AIRY_H_ */
