#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "integration.h"

int
gsl_integration_qag (const gsl_function *f,
		     double a, double b,
		     double epsabs, double epsrel, size_t limit,
		     int key,
		     gsl_integration_workspace * workspace,
		     double * result, double * abserr)
{
  int status ;
  gsl_integration_rule * integration_rule = gsl_integration_qk15 ;

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
      integration_rule = gsl_integration_qk15 ;
      break ;
    case GSL_INTEG_GAUSS21:
      integration_rule = gsl_integration_qk21 ;
      break ;
    case GSL_INTEG_GAUSS31:
      integration_rule = gsl_integration_qk31 ; 
      break ;
    case GSL_INTEG_GAUSS41:
      integration_rule = gsl_integration_qk41 ;
      break ;      
    case GSL_INTEG_GAUSS51:
      integration_rule = gsl_integration_qk51 ;
      break ;      
    case GSL_INTEG_GAUSS61:
      integration_rule = gsl_integration_qk61 ;
      break ;      
    default:
      GSL_ERROR("value of key does specify a known integration rule", 
		GSL_EINVAL) ;
    }

  status = gsl_integration_qag_impl (f, a, b, epsabs, epsrel, limit,
				     workspace, 
				     result, abserr, 
				     integration_rule) ;

  return status ;
}

