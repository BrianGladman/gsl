#include <stdio.h>
#include <math.h>
#include <gsl_integration.h>

double f (double x) ;

int main (void)
{
  double result, abserr, resabs, resasc ;

  gsl_integration_qk15(f,0.0,1.0,&result,&abserr,&resabs,&resasc) ;

  printf("15: result = %.30g, abserr = %.18g, resabs = %.18g, resasc = %.18g\n",
	 result, abserr, resabs, resasc) ;

  gsl_integration_qk21(f,0.0,1.0,&result,&abserr,&resabs,&resasc) ;

  printf("21: result = %.18g, abserr = %.18g, resabs = %.18g, resasc = %.18g\n",
	 result, abserr, resabs, resasc) ;

  gsl_integration_qk31(f,0.0,1.0,&result,&abserr,&resabs,&resasc) ;

  printf("31: result = %.18g, abserr = %.18g, resabs = %.18g, resasc = %.18g\n",
	 result, abserr, resabs, resasc) ;

  gsl_integration_qk41(f,0.0,1.0,&result,&abserr,&resabs,&resasc) ;

  printf("41: result = %.18g, abserr = %.18g, resabs = %.18g, resasc = %.18g\n",
	 result, abserr, resabs, resasc) ;

  gsl_integration_qk51(f,0.0,1.0,&result,&abserr,&resabs,&resasc) ;

  printf("51: result = %.18g, abserr = %.18g, resabs = %.18g, resasc = %.18g\n",
	 result, abserr, resabs, resasc) ;


  gsl_integration_qk61(f,0.0,1.0,&result,&abserr,&resabs,&resasc) ;

  printf("61: result = %.18g, abserr = %.18g, resabs = %.18g, resasc = %.18g\n",
	 result, abserr, resabs, resasc) ;



  return 0 ;
} 

double f (double x) {
  return  x*x ;
}
