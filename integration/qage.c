#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_integration.h>

int
gsl_integration_qage (double (*f)(double x),
		      double a, double b,
		      double epsabs, double epsrel,
		      int key,
		      gsl_integration_workspace * workspace,
		      size_t * last,
		      double * result, double * abserr, size_t * neval)
{
  int status ;
  size_t nqeval = 0;
  gsl_integration_rule_t * integration_rule = &gsl_integration_qk15 ;

  if (key < GSL_INTEG_GAUSS15)
    {
      key = GSL_INTEG_GAUSS15 ;
    } 
  else if (key > GSL_INTEG_GAUSS61) 
    {
      key = GSL_INTEG_GAUSS61 ;
    }

  switch (key) 
    {
    case GSL_INTEG_GAUSS15:
      integration_rule = &gsl_integration_qk15 ;
      break ;
    case GSL_INTEG_GAUSS21:
      integration_rule = &gsl_integration_qk21 ;
      break ;
    case GSL_INTEG_GAUSS31:
      integration_rule = &gsl_integration_qk31 ; 
      break ;
    case GSL_INTEG_GAUSS41:
      integration_rule = &gsl_integration_qk41 ;
      break ;      
    case GSL_INTEG_GAUSS51:
      integration_rule = &gsl_integration_qk51 ;
      break ;      
    case GSL_INTEG_GAUSS61:
      integration_rule = &gsl_integration_qk61 ;
      break ;      
    default:
      GSL_ERROR("value of key does specify a known integration rule", 
		GSL_EINVAL) ;
    }

  status = gsl_integration_qage_impl (f, a, b, epsabs, epsrel, 
				      workspace, last, 
				      result, abserr, &nqeval, 
				      integration_rule) ;

  /* convert from number of quadrature rule evaluations to number of
     function evaluations */

  if (key == GSL_INTEG_GAUSS15)
    {
      *neval = 15 * nqeval ;
    }
  else 
    {
      *neval = (10 * key + 1) * nqeval ;
    }

  return status ;
}

