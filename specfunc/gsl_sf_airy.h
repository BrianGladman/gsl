/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_AIRY_H__
#define __GSL_SF_AIRY_H__

#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf_result.h>

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


/* Airy function Ai(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_airy_Ai_impl(const double x, const gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Ai_e(const double x, gsl_mode_t mode, gsl_sf_result * result);


/* Airy function Bi(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_airy_Bi_impl(const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Bi_e(const double x, gsl_mode_t mode, gsl_sf_result * result);


/* scaled Ai(x):
 *                     Ai(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Ai_scaled_impl(const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Ai_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);


/* scaled Bi(x):
 *                     Bi(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_scaled_impl(const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Bi_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);


/* derivative Ai'(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_airy_Ai_deriv_impl(const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Ai_deriv_e(const double x, gsl_mode_t mode, gsl_sf_result * result);


/* derivative Bi'(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_airy_Bi_deriv_impl(const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Bi_deriv_e(const double x, gsl_mode_t mode, gsl_sf_result * result);


/* scaled derivative Ai'(x):
 *                     Ai'(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai'(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Ai_deriv_scaled_impl(const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Ai_deriv_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);


/* scaled derivative:
 *                     Bi'(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi'(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_deriv_scaled_impl(const double x, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_airy_Bi_deriv_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);


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


__END_DECLS

#endif /* __GSL_SF_AIRY_H__ */
