#include <gsl_integration.h>
#include "qk31.h"

int
gsl_integration_qk31 (double (*f) (double x),
		      double a, double b,
		      double * result, double * abserr,
		      double * resabs, double * resasc)
{
  double fv1[16], fv2[16];
  gsl_integration_qk(16,xgk,wg,wgk,fv1,fv2,
		     f,a,b,result,abserr,resabs,resasc) ;
  return 0 ;
}
