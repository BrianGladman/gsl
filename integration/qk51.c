#include <config.h>
#include <gsl_integration.h>
#include "qk51.h"

void
gsl_integration_qk51 (const gsl_function *f,
		      double a, double b,
		      double * result, double * abserr,
		      double * resabs, double * resasc)
{
  double fv1[26], fv2[26];
  gsl_integration_qk(26,xgk,wg,wgk,fv1,fv2,
		     f,a,b,result,abserr,resabs,resasc) ;
}
