#include <gsl_integration.h>
#include "qk21.h"

int
gsl_integration_qk21 (double (*f) (double x),
		      double a, double b,
		      double * result, double * abserr,
		      double * resabs, double * resasc)
{
  double fv1[11], fv2[11];
  gsl_integration_qk(11,xgk,wg,wgk,fv1,fv2,
		     f,a,b,result,abserr,resabs,resasc) ;
  return 0 ;
}
