#include <stdio.h>
#include <math.h>
#include <gsl_integration.h>

double qtest1 (double x) ;
double g (double x) ;
double stepfn (double x) ;

double alpha = 0 ;

int main (void)
{
  size_t neval ;
  int status ;
      
  {
    double result = 0, abserr = 0, resabs = 0, resasc = 0 ;
    quad_rule (f,0.0,1.0,&result,&abserr,&resabs,&resasc) ;
  }

  gsl_integration_qng(f,0.0,1.0,0.0,1e-5,&result,&abserr,&neval) ;
  
  printf("qng: result = %.18g, abserr = %.18g, neval = %d\n",
	 result, abserr, neval) ;

#ifdef JUNK
  {
    double alist[1000], blist[1000], rlist[1000], elist[1000];
    size_t iord[1000] ;
    size_t last;
    result = 0 ; abserr=0; neval=0  ;
    gsl_integration_qage(f, 0.0, 1.0, 0, 1e-10, 6, 1000,
			 alist, blist, rlist, elist, iord, &last,
			 &result, &abserr, &neval) ;
    printf("qage: result = %.18g, abserr = %.18g, neval = %d\n",
	   result, abserr, neval) ;
  }
#endif

  {
    double alist[1000], blist[1000], rlist[1000], elist[1000];
    size_t iord[1000] ;
    size_t last;
    result = 0 ; abserr=0; neval=0  ;
    status = gsl_integration_qagse(f, 0.0, 1.0, 1e-20, 1e-17, 10,
			 alist, blist, rlist, elist, iord, &last,
			 &result, &abserr, &neval) ;
    printf("qagse: result = %.18g, abserr = %.18g, neval = %d\n",
	   result, abserr, neval) ;
    printf("status=%d\n",status) ;
  }

#ifdef JUNK
  {
    double alist[1000], blist[1000], rlist[1000], elist[1000];
    size_t iord[1000] ;
    size_t last;
    result = 0 ; abserr=0; neval=0  ;
    status = gsl_integration_qagse(g, 0.0, 1.0, 0, 1e-10, 10,
			 alist, blist, rlist, elist, iord, &last,
			 &result, &abserr, &neval) ;
    printf("qagse: result = %.18g, abserr = %.18g, neval = %d\n",
	   result, abserr, neval) ;
    printf("status=%d\n",status) ;
  }


  {
    double alist[1000], blist[1000], rlist[1000], elist[1000];
    size_t iord[1000] ;
    size_t last;
    result = 0 ; abserr=0; neval=0  ;
    status = gsl_integration_qagse(stepfn, 0.0, 1.0, 0, 1e-10, 3,
				   alist, blist, rlist, elist, iord, &last,
				   &result, &abserr, &neval) ;
    printf("qagse: result = %.18g, abserr = %.18g, neval = %d\n",
	   result, abserr, neval) ;
    printf("status=%d\n",status) ;
  }
#endif

  return 0 ;
} 

double qtest1 (double x) {
  return pow(x,alpha) * log(1/x) ;
}

double g (double x) {
  return 1/sqrt(fabs(x*x + 2*x - 2)) ;
}


double stepfn (double x) {
  if (x < 1.0/3.0)
    {
      return 0 ;
    }
  
  return 1 ;
}




