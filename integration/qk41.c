#include <gsl_integration.h>
#include "qk41.h"

void
gsl_integration_qk41 (double (*f) (double x),
		      double a, double b,
		      double * result, double * abserr,
		      double * resabs, double * resasc)
{
  double fv1[21], fv2[21];
  gsl_integration_qk(21,xgk,wg,wgk,fv1,fv2,
		     f,a,b,result,abserr,resabs,resasc) ;
}
