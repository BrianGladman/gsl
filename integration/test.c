#include <stdio.h>
#include <math.h>
#include <gsl_integration.h>

double f (double x) ;

int main (void)
{
  double result, abserr, resabs, resasc ;
  int neval ;

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


  gsl_integration_qng(f,0.0,1.0,0.0,1e-5,&result,&abserr,&neval) ;

  printf("qnd: result = %.18g, abserr = %.18g, neval = %d\n",
	 result, abserr, neval) ;



  return 0 ;
} 

double f (double x) {
  return  sin(60*M_PI*x) ;
}
