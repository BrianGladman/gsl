#include <config.h>
#include <gsl_integration.h>
#include "qk15.h"

void
gsl_integration_qk15 (double (*f) (double x),
		      double a, double b,
		      double * result, double * abserr,
		      double * resabs, double * resasc)
{
  double fv1[8], fv2[8];
  gsl_integration_qk(8,xgk,wg,wgk,fv1,fv2,
		     f,a,b,result,abserr,resabs,resasc) ;
}
