/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef AIRY_IMPL_H_
#define AIRY_IMPL_H_


int gsl_sf_airy_Ai_impl(double x, double * result);
int gsl_sf_airy_Bi_impl(double x, double * result);
int gsl_sf_airy_Ai_scaled_impl(double x, double * result);
int gsl_sf_airy_Bi_scaled_impl(double x, double * result);

int gsl_sf_airy_Ai_deriv_scaled_impl(double x, double * result);
int gsl_sf_airy_Bi_deriv_scaled_impl(double x, double * result);
int gsl_sf_airy_Ai_deriv_impl(double x, double * result);
int gsl_sf_airy_Bi_deriv_impl(double x, double * result);


#endif /* !AIRY_IMPL_H_ */
