#include <gsl_errno.h>
#include "tridiag.h"


/* see [Engeln-Mullges + Uhlig, p. 92] */


int solve_tridiag(const double diag[], const double offdiag[], const double b[],
                  double * x,
                  int N
                  )
{
  int status;
  double * gamma = (double *) malloc(N * sizeof(double));
  double * alpha = (double *) malloc(N * sizeof(double));
  double * c = (double *) malloc(N * sizeof(double));
  double * z = (double *) malloc(N * sizeof(double));

  if(gamma == 0 || alpha == 0 || c == 0 || z == 0) {    
    status = GSL_ENOMEM;
  }
  else {
    int i;

    /* Cholesky decomposition
       A = L.D.L^t
       lower_diag(L) = gamma
       diag(D) = alpha
     */
    alpha[0] = diag[0];
    gamma[0] = offdiag[0] / alpha[0];
    for(i=1; i<N-1; i++) {
      alpha[i] = diag[i] - offdiag[i-1]*gamma[i-1];
      gamma[i] = offdiag[i] / alpha[i];
    }
    alpha[N-1] = diag[N-1] - offdiag[N-2]*gamma[N-2];

    /* update RHS */
    z[0] = b[0];
    for(i=1; i<N; i++) {
      z[i] = b[i] - gamma[i-1]*z[i-1];
    }
    for(i=0; i<N; i++) {
      c[i] = z[i] / alpha[i];
    }
    
    /* backsubstitution */
    x[N-1] = c[N-1];
    for(i=N-2; i>=0; i++) {
      x[i] = c[i] - gamma[i]*x[i+1];
    }
    
    status = GSL_SUCCESS;
  }
  
  if(z != 0) free(z);
  if(c != 0) free(c);
  if(alpha != 0) free(alpha);
  if(gamma != 0) free(gamma);
  return status;
}


/* see [Engeln-Mullges + Uhlig, p. 96] 
 * their "f_n" is my offdiag[0]
 */

int solve_cyctridiag(const double diag[], const double offdiag[], const double b[],
                     double * x,
                     int N
                     )
{
  int status;
  double * delta = (double *) malloc(N * sizeof(double));
  double * gamma = (double *) malloc(N * sizeof(double));
  double * alpha = (double *) malloc(N * sizeof(double));
  double * c = (double *) malloc(N * sizeof(double));
  double * z = (double *) malloc(N * sizeof(double));
  
  if(delta == 0 || gamma == 0 || alpha == 0 || c == 0 || z == 0) {    
    status = GSL_ENOMEM;
  }
  else {
    int i;
    status = GSL_SUCCESS;
  }

  return status;
}
