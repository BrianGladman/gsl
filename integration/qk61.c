#include <gsl_integration.h>
#include "qk61.h"

int
gsl_integration_qk61 (double (*f) (double x),
		      double a, double b,
		      double * result, double * abserr,
		      double * resabs, double * resasc)
{
  double fv1[31], fv2[31];
  gsl_integration_qk(31,xgk,wg,wgk,fv1,fv2,
		     f,a,b,result,abserr,resabs,resasc) ;
  return 0 ;
}
