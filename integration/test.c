#include <stdio.h>
#include <gsl_integration.h>

double f (double x) ;

int main (void)
{
  double result, abserr, resabs, resasc ;

  gsl_integration_qk15(f,0,0.001,&result,&abserr,&resabs,&resasc) ;

  printf("result = %.18g, abserr = %.18g, resabs = %.18g, resasc = %.18g\n",
	 result, abserr, resabs, resasc) ;

} ;

double f (double x) {
  return x*x ;
}
