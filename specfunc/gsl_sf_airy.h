/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_AIRY_H__
#define __GSL_SF_AIRY_H__

#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf_result.h>


/* Airy function Ai(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_airy_Ai_impl(double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Ai_e(double x, gsl_mode_t mode, gsl_sf_result * result);


/* Airy function Bi(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_airy_Bi_impl(double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Bi_e(double x, gsl_mode_t mode, gsl_sf_result * result);


/* scaled Ai(x):
 *                     Ai(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Ai_scaled_impl(double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Ai_scaled_e(double x, gsl_mode_t mode, gsl_sf_result * result);


/* scaled Bi(x):
 *                     Bi(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_scaled_impl(double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Bi_scaled_e(double x, gsl_mode_t mode, gsl_sf_result * result);


/* derivative Ai'(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_airy_Ai_deriv_impl(double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Ai_deriv_e(double x, gsl_mode_t mode, gsl_sf_result * result);


/* derivative Bi'(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_airy_Bi_deriv_impl(double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Bi_deriv_e(double x, gsl_mode_t mode, gsl_sf_result * result);


/* scaled derivative Ai'(x):
 *                     Ai'(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai'(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Ai_deriv_scaled_impl(double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Ai_deriv_scaled_e(double x, gsl_mode_t mode, gsl_sf_result * result);


/* scaled derivative:
 *                     Bi'(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi'(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_deriv_scaled_impl(double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Bi_deriv_scaled_e(double x, gsl_mode_t mode, gsl_sf_result * result);


/* Zeros of Ai(x)
 */
int gsl_sf_airy_zero_Ai_impl(int s, gsl_sf_result * result);
int gsl_sf_airy_zero_Ai_e(int s, gsl_sf_result * result);


/* Zeros of Bi(x)
 */
int gsl_sf_airy_zero_Bi_impl(int s, gsl_sf_result * result);
int gsl_sf_airy_zero_Bi_e(int s, gsl_sf_result * result);


/* Zeros of Ai'(x)
 */
int gsl_sf_airy_zero_Ai_deriv_impl(int s, gsl_sf_result * result);
int gsl_sf_airy_zero_Ai_deriv_e(int s, gsl_sf_result * result);


/* Zeros of Bi'(x)
 */
int gsl_sf_airy_zero_Bi_deriv_impl(int s, gsl_sf_result * result);
int gsl_sf_airy_zero_Bi_deriv_e(int s, gsl_sf_result * result);



#endif /* __GSL_SF_AIRY_H__ */
