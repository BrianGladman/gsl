#include <gsl_multimin.h>
#include <gsl_blas_types.h>
#include <gsl_blas.h>

int
gsl_multimin_test_gradient_sqr_norm(gsl_multimin_fdf_history *h,double epsabs)
{
  double sqr_norm;

  if (epsabs < 0.0)
    GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);

  gsl_blas_ddot(h->g,h->g,&sqr_norm);
  
  if (sqr_norm<epsabs)
    return GSL_SUCCESS;

  return GSL_CONTINUE;
}
