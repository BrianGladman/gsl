#include <config.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <gsl_math.h>
#include <gsl_integration.h>

void
gsl_integration_qelg (size_t * n, double epstab[], 
		      double * result, double * abserr,
		      double res3la[], size_t * nres)
{
  const double current = epstab[*n] ;

  double absolute = GSL_DBL_MAX ;
  double relative = 5 * GSL_DBL_EPSILON * fabs(current) ;

  const size_t newelm = (*n)/2 ;
  const size_t n_orig = (*n) ;
  size_t n_final = (*n) ;
  size_t i ;
  const size_t nres_orig = *nres ;
  *result = current ;
  *nres = (*nres) + 1 ;
  *abserr = GSL_DBL_MAX ;

#ifdef DEBUG
  for (i= 0; i<32; i++) {
    printf("QELG: TAB i = %d epstab(i) %.18e\n",i+1,epstab[i])  ;
  } ;
#endif

  if ((*n) < 2) 
    {
      *result = current ;
      *abserr = GSL_MAX_DBL (absolute,relative) ;
      return ;
    }
  
  epstab[(*n)+2] = epstab[(*n)] ;
  epstab[(*n)] = GSL_DBL_MAX ;

#ifdef DEBUG  
  printf("QELG: outside loop, newelm = %d, n_orig = %d\n", newelm, n_orig) ;
#endif

  for (i = 0; i < newelm ; i++)
    {
      double res = epstab[(*n) - 2*i + 2] ;
      double e0 = epstab[(*n) - 2*i - 2] ;
      double e1 = epstab[(*n) - 2*i - 1] ;
      double e2 = res ;
      
      double e1abs = fabs(e1) ;
      double delta2 = e2 - e1 ;
      double err2 = fabs(delta2) ;
      double tol2 = GSL_MAX_DBL(fabs(e2),e1abs)*GSL_DBL_EPSILON ;
      double delta3 = e1 - e0 ;
      double err3 = fabs(delta3) ;
      double tol3 = GSL_MAX_DBL(e1abs,fabs(e0))*GSL_DBL_EPSILON ;
      
      double e3, delta1, err1, tol1, ss ;

#ifdef DEBUG  
      printf("QELG: in loop, i = %d, newelm = %d\n", i,newelm) ;
#endif

      if (err2 <= tol2 && err3 <= tol3)
        {
          /* If e0, e1 and e2 are equal to within machine accuracy,
           convergence is assumed.  */
#ifdef DEBUG  
	  printf("QELG: err2 <= tol2 && err3 <= tol3\n") ;
#endif
          *result = res ;
          absolute = err2 + err3 ;
          relative = 5 * GSL_DBL_EPSILON * fabs(res) ;
	  *abserr = GSL_MAX_DBL(absolute, relative) ;
          return ;
        }

      e3 = epstab[(*n) - 2*i] ;
      epstab[(*n) - 2*i] = e1 ;
      delta1 = e1 - e3 ;
      err1 = fabs(delta1) ;
      tol1 = GSL_MAX_DBL(e1abs, fabs(e3)) * GSL_DBL_EPSILON ;
      
      /* If two elements are very close to each other, omit a part of
         the table by adjusting the value of n */
#ifdef DEBUG        
      printf("QELG: err1 = %.18e tol1 = %.18e\n",err1,tol1) ;
      printf("QELG: err2 = %.18e tol2 = %.18e\n",err2,tol2) ;
      printf("QELG: err3 = %.18e tol3 = %.18e\n",err3,tol3) ;
#endif
      if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3)
        {
#ifdef DEBUG  	  
	  printf("QELG: err1 <= tol1 || err2 <= tol2 || err3 <= tol3\n") ;
#endif
          n_final = 2*i ;
          break ;
        }
      
      ss = (1/delta1 + 1/delta2) - 1/delta3 ;

#ifdef DEBUG  
      printf("QELG: ss = %.18e\n",ss) ;
#endif
      /* Test to detect irregular behaviour in the table, and
         eventually omit a part of the table by adjusting the value of
         n. */

      if (fabs(ss*e1) <= 0.0001) 
        {
#ifdef DEBUG	 
	  printf("QELG: fabs(ss*e1) <= 0.0001\n");
#endif
          n_final = 2*i ;
          break ;
        }

      /* Compute a new element and eventually adjust the value of
         result. */
      
      res = e1 + 1/ss ;
      epstab[(*n) - 2*i] = res ;

      {
	const double error = err2 + fabs(res - e2) + err3 ;
#ifdef DEBUG  
	printf("QELG: *abserr = %.18e, error = %.18e\n", *abserr, error) ;
#endif
	if (error <= *abserr) 
	  {
#ifdef DEBUG	    
	    printf("QELG: setting *abserr to %.18e\n",error) ;
#endif	    
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
#ifdef DEBUG  
	printf("QELG: nfinal reduced to %d\n",n_final) ;
#endif
      }
  }
  
  if (n_orig % 2 == 1)
    {
      for (i = 0 ; i <= newelm ; i++)
	{
#ifdef DEBUG	  
	  printf("QELG: A:copying epstab[%d] into epstab[%d]\n",i*2+3,i*2+1) ;
#endif	  
	  epstab[1+i*2] = epstab[i*2+3] ;
	}
    }
  else
    {
      for (i = 0 ; i <= newelm ; i++)
	{
#ifdef DEBUG	  
	  printf("QELG: A:copying epstab[%d] into epstab[%d]\n",i*2+2,i*2) ;
#endif
	  epstab[i*2] = epstab[i*2+2] ;
	}
    }

  if (n_orig != n_final) {
    printf ("n_orig = %d, n_final = %d\n",n_orig,n_final) ;
    for (i = 0 ; i <= n_final ; i++)
      {
#ifdef DEBUG  
	printf("QELG: B:copying epstab[%d] into epstab[%d]\n",n_orig-n_final+i,i) ;
#endif
        epstab[i] = epstab[n_orig - n_final + i] ;
      }
  }

  *n = n_final ;

  if (nres_orig < 3) 
    {
#ifdef DEBUG  
      printf("QELG: setting term %d to %.18e\n",nres_orig, *result) ;
#endif
      res3la[nres_orig] = *result ;
      *abserr = GSL_DBL_MAX ;
    } 
  else 
    {  /* Compute error estimate */
      *abserr = (fabs(*result - res3la[2]) + fabs(*result - res3la[1])
		+ fabs(*result - res3la[0])) ;
#ifdef DEBUG  
      printf("QELG: error estimate computed as %.18e\n",*abserr) ;
      printf("QELG: result = %.18e\n", *result) ;
      printf("QELG: term1 = %.18e\n", res3la[2]) ;
      printf("QELG: term2 = %.18e\n", res3la[1]) ;
      printf("QELG: term3 = %.18e\n", res3la[0]) ;
#endif
      res3la[0] = res3la[1] ;
      res3la[1] = res3la[2] ;
      res3la[2] = *result ;
    }

  *abserr = GSL_MAX_DBL(*abserr, 5 * GSL_DBL_EPSILON * fabs(*result)) ;

#ifdef DEBUG  
      printf("QELG: abserr = %.18e\n",*abserr) ;
#endif

  return ;
}
    


