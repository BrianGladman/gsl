/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_LOG_H_
#define GSL_SF_LOG_H_


/* Provide a logarithm function with GSL semantics.
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_log_impl(double x, double * result);
int     gsl_sf_log_e(double x, double * result);
double  gsl_sf_log(double x);


/* Log(|x|)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_log_abs_impl(double x, double * result);
int     gsl_sf_log_abs_e(double x, double * result);
double  gsl_sf_log_abs(double x);


/* Complex Logarithm
 *   exp(lnr + I theta) = zr + I zi
 * Returns argument in [-pi,pi].
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_log_impl(double zr, double zi, double * lnr, double * theta);
int gsl_sf_complex_log_e(double zr, double zi, double * lnr, double * theta);


/* Log(1 + x)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_log_1plusx_impl(double x, double * result);
int     gsl_sf_log_1plusx_e(double x, double * result);
double  gsl_sf_log_1plusx(double x);


#ifdef HAVE_INLINE
extern inline
int
gsl_sf_log_impl(const double x, double * result)
{
  if(x <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    *result = log(x);
    return GSL_SUCCESS;
  }
}
extern inline
int
gsl_sf_log_abs_impl(const double x, double * result)
{
  if(x == 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    *result = log(fabs(x));
    return GSL_SUCCESS;
  }
}
#endif  /* HAVE_INLINE */


#endif /* !GSL_SF_LOG_H_ */
