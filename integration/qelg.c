#include <stdio.h>

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <gsl_integration.h>

#include "max.h"

void
gsl_integration_qelg (size_t * n, double epstab[], 
		      double * result, double * abserr,
		      double res3la[], size_t * nres)
{
  const double current = epstab[*n] ;

  double absolute = DBL_MAX ;
  double relative = 5 * DBL_EPSILON * fabs(current) ;

  const size_t newelm = (*n)/2 ;
  const size_t n_orig = (*n) ;
  size_t n_final = (*n) ;
  size_t i ;
  const size_t nres_orig = *nres ;
  *result = current ;
  *nres = (*nres) + 1 ;
  *abserr = DBL_MAX ;
  if ((*n) < 2) 
    {
      *result = current ;
      *abserr = max(absolute,relative) ;
      return ;
    }
  
  epstab[(*n)+2] = epstab[(*n)] ;
  epstab[(*n)] = DBL_MAX ;
  
  printf("QELG: outside loop, newelm = %d, n_orig = %d\n", newelm, n_orig) ;

  for (i = 0; i < newelm ; i++)
    {
      double res = epstab[(*n) - 2*i + 2] ;
      double e0 = epstab[(*n) - 2*i - 2] ;
      double e1 = epstab[(*n) - 2*i - 1] ;
      double e2 = res ;
      
      double e1abs = fabs(e1) ;
      double delta2 = e2 - e1 ;
      double err2 = fabs(delta2) ;
      double tol2 = max(fabs(e2),e1abs)*DBL_EPSILON ;
      double delta3 = e1 - e0 ;
      double err3 = fabs(delta3) ;
      double tol3 = max(e1abs,fabs(e0))*DBL_EPSILON ;
      
      double e3, delta1, err1, tol1, ss ;

      printf("QELG: in loop, i = %d, newelm = %d\n", i,newelm) ;

      if (err2 <= tol2 && err3 <= tol3)
        {
          /* If e0, e1 and e2 are equal to within machine accuracy,
           convergence is assumed.  */

	  printf("QELG: err2 <= tol2 && err3 <= tol3\n") ;

          *result = res ;
          absolute = err2 + err3 ;
          relative = 5 * DBL_EPSILON * fabs(res) ;
	  *abserr = max(absolute, relative) ;
          return ;
        }

      e3 = epstab[(*n) - 2*i] ;
      epstab[(*n) - 2*i] = e1 ;
      delta1 = e1 - e3 ;
      err1 = fabs(delta1) ;
      tol1 = max(e1abs, fabs(e3)) * DBL_EPSILON ;
      
      /* If two elements are very close to each other, omit a part of
         the table by adjusting the value of n */
      
      printf("QELG: err1 = %g tol1 = %g\n",err1,tol1) ;
      printf("QELG: err2 = %g tol2 = %g\n",err2,tol2) ;
      printf("QELG: err3 = %g tol3 = %g\n",err3,tol3) ;

      if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3)
        {
	  printf("QELG: err1 <= tol1 || err2 <= tol2 || err3 <= tol3\n") ;
          n_final = 2*i ;
          break ;
        }
      
      ss = (1/delta1 + 1/delta2) - 1/delta3 ;

      printf("QELG: ss = %g\n",ss) ;

      /* Test to detect irregular behaviour in the table, and
         eventually omit a part of the table by adjusting the value of
         n. */

      if (fabs(ss*e1) <= 0.0001) 
        {
	  printf("QELG: fabs(ss*e1) <= 0.0001\n");
          n_final = 2*i ;
          break ;
        }

      /* Compute a new element and eventually adjust the value of
         result. */
      
      res = e1 + 1/ss ;
      epstab[(*n) - 2*i] = res ;

      {
	const double error = err2 + fabs(res - e2) + err3 ;
	printf("QELG: *abserr = %g, error = %g\n", *abserr, error) ;

	if (error <= *abserr) 
	  {
	    printf("QELG: setting *abserr to %g\n",error) ;
	    *abserr = error ;
	    *result = res ;
	  } 
      }
    }

  /* Shift the table */
  
  {
    const size_t limexp = 50 - 1 ;

    if (n_final == limexp)
      {
	n_final = 2 * (limexp/2)  ;
	printf("QELG: nfinal reduced to %d\n",n_final) ;
      }
  }
  
  if (n_orig % 2 == 1)
    {
      for (i = 0 ; i <= newelm ; i++)
	{
	  printf("QELG: A:copying epstab[%d] into epstab[%d]\n",i*2+3,i*2+1) ;
	  epstab[1+i*2] = epstab[i*2+3] ;
	}
    }
  else
    {
      for (i = 0 ; i <= newelm ; i++)
	{
	  printf("QELG: A:copying epstab[%d] into epstab[%d]\n",i*2+2,i*2) ;
	  epstab[i*2] = epstab[i*2+2] ;
	}
    }

  if (n_orig != n_final) {
    printf ("n_orig = %d, n_final = %d\n",n_orig,n_final) ;
    for (i = 0 ; i <= n_final ; i++)
      {
	printf("QELG: B:copying epstab[%d] into epstab[%d]\n",n_orig-n_final+i,i) ;
        epstab[i] = epstab[n_orig - n_final + i] ;
      }
  }

  *n = n_final ;

  if (nres_orig < 3) 
    {
      printf("QELG: setting term %d to %g\n",nres_orig, *result) ;
      res3la[nres_orig] = *result ;
      *abserr = DBL_MAX ;
    } 
  else 
    {  /* Compute error estimate */
      *abserr = (fabs(*result - res3la[2]) + fabs(*result - res3la[1])
		+ fabs(*result - res3la[0])) ;
      printf("QELG: error estimate computed as %g\n",*abserr) ;
      printf("QELG: result = %g\n", *result) ;
      printf("QELG: term1 = %g\n", res3la[2]) ;
      printf("QELG: term2 = %g\n", res3la[1]) ;
      printf("QELG: term3 = %g\n", res3la[0]) ;
      res3la[0] = res3la[1] ;
      res3la[1] = res3la[2] ;
      res3la[2] = *result ;
    }

  *abserr = max(*abserr, 5 * DBL_EPSILON * fabs(*result)) ;

#ifdef JUN
  for (i= 0; i<32; i++) {
    printf("QELG: TAB i = %d, epstab[i] = %g\n",i,epstab[i])  ;
  } ;
#endif
  return ;
}
    


