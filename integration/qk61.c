#include <config.h>
#include <gsl_integration.h>
#include "qk61.h"

void
gsl_integration_qk61 (const gsl_function *f,
		      double a, double b,
		      double * result, double * abserr,
		      double * resabs, double * resasc)
{
  double fv1[31], fv2[31];
  gsl_integration_qk(31,xgk,wg,wgk,fv1,fv2,
		     f,a,b,result,abserr,resabs,resasc) ;
}
