


/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdlib.h>
#include <gsl_errno.h>
#include "tridiag.h"


/* see [Engeln-Mullges + Uhlig, p. 92] */


int 
solve_tridiag (const double diag[], const double offdiag[], const double b[],
	       double *x,
	       size_t N)
{
  int status;
  double *gamma = (double *) malloc (N * sizeof (double));
  double *alpha = (double *) malloc (N * sizeof (double));
  double *c = (double *) malloc (N * sizeof (double));
  double *z = (double *) malloc (N * sizeof (double));

  if (gamma == 0 || alpha == 0 || c == 0 || z == 0)
    {
      status = GSL_ENOMEM;
    }
  else
    {
      size_t i, j;

      /* Cholesky decomposition
         A = L.D.L^t
         lower_diag(L) = gamma
         diag(D) = alpha
       */
      alpha[0] = diag[0];
      gamma[0] = offdiag[0] / alpha[0];
      for (i = 1; i < N - 1; i++)
	{
	  alpha[i] = diag[i] - offdiag[i - 1] * gamma[i - 1];
	  gamma[i] = offdiag[i] / alpha[i];
	}
      alpha[N - 1] = diag[N - 1] - offdiag[N - 2] * gamma[N - 2];

      /* update RHS */
      z[0] = b[0];
      for (i = 1; i < N; i++)
	{
	  z[i] = b[i] - gamma[i - 1] * z[i - 1];
	}
      for (i = 0; i < N; i++)
	{
	  c[i] = z[i] / alpha[i];
	}

      /* backsubstitution */
      x[N - 1] = c[N - 1];
      if (N >= 2)
	{
	  for (i = N - 2, j = 0; j <= N - 2; j++, i--)
	    {
	      x[i] = c[i] - gamma[i] * x[i + 1];
	    }
	}

      status = GSL_SUCCESS;
    }

  if (z != 0)
    free (z);
  if (c != 0)
    free (c);
  if (alpha != 0)
    free (alpha);
  if (gamma != 0)
    free (gamma);
  return status;
}


/* see [Engeln-Mullges + Uhlig, p. 96] 
 * their "f_n" is my offdiag[0]
 */

int 
solve_cyctridiag (const double diag[], const double offdiag[], const double b[],
		  double *x,
		  size_t N
)
{
  int status;
  double *delta = (double *) malloc (N * sizeof (double));
  double *gamma = (double *) malloc (N * sizeof (double));
  double *alpha = (double *) malloc (N * sizeof (double));
  double *c = (double *) malloc (N * sizeof (double));
  double *z = (double *) malloc (N * sizeof (double));

  if (delta == 0 || gamma == 0 || alpha == 0 || c == 0 || z == 0)
    {
      status = GSL_ENOMEM;
    }
  else
    {
      size_t i, j;
      double sum = 0.0;

      /* factor */
      alpha[0] = diag[0];
      gamma[0] = offdiag[1] / alpha[0];
      delta[0] = offdiag[0] / alpha[0];
      for (i = 1; i < N - 2; i++)
	{
	  alpha[i] = diag[i] - offdiag[i] * gamma[i - 1];
	  gamma[i] = offdiag[i + 1] / alpha[i];
	  delta[i] = -delta[i - 1] * offdiag[i] / alpha[i];
	}
      for (i = 0; i < N - 2; i++)
	{
	  sum += alpha[i] * delta[i] * delta[i];
	}
      alpha[N - 2] = diag[N - 2] - offdiag[N - 2] * gamma[N - 3];
      gamma[N - 2] = (offdiag[N - 1] - offdiag[N - 2] * delta[N - 3]) / alpha[N - 2];
      alpha[N - 1] = diag[N - 1] - sum - offdiag[N - 1] * gamma[N - 2] * gamma[N - 2];

      /* update */
      z[0] = b[0];
      for (i = 1; i < N - 1; i++)
	{
	  z[i] = b[i] - z[i - 1] * gamma[i - 1];
	}
      sum = 0.0;
      for (i = 0; i < N - 2; i++)
	{
	  sum += delta[i] * z[i];
	}
      z[N - 1] = b[N - 1] - sum - gamma[N - 2] * z[N - 2];
      for (i = 0; i < N; i++)
	{
	  c[i] = z[i] / alpha[i];
	}

      /* backsubstitution */
      x[N - 1] = c[N - 1];
      x[N - 2] = c[N - 2] - gamma[N - 2] * x[N - 1];
      if (N >= 3)
	{
	  for (i = N - 3, j = 0; j <= N - 3; j++, i--)
	    {
	      x[i] = c[i] - gamma[i] * x[i + 1] - delta[i] * x[N - 1];
	    }
	}

      status = GSL_SUCCESS;
    }

  return status;
}
