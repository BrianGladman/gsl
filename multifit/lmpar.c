/* multifit/lmpar.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <gsl/gsl_permute_vector_double.h>

#include "qrsolv.c"

static size_t
count_nsing (const gsl_matrix * r)
{
  /* Count the number of nonsingular entries. Returns the index of the
     first entry which is singular. */

  size_t n = r->size2;
  size_t i;

  for (i = 0; i < n; i++)
    {
      double rii = gsl_matrix_get (r, i, i);

      if (rii == 0)
	{
	  break;
	}
    }

  return i;
}


static void
compute_newton_direction (const gsl_matrix * r, const gsl_permutation * perm,
			  const gsl_vector * qtf, gsl_vector * x)
{

  /* Compute and store in x the Gauss-Newton direction. If the
     Jacobian is rank-deficient then obtain a least squares
     solution. */

  const size_t n = r->size2;
  size_t i, j, nsing;

  gsl_vector_memcpy (x, qtf);

  nsing = count_nsing (r);

  for (i = nsing; i < n; i++)
    {
      gsl_vector_set (x, i, 0.0);
    }

  for (j = nsing - 1; j > 0; j--)
    {
      double rjj = gsl_matrix_get (r, j, j);
      double temp = gsl_vector_get (x, j) / rjj;

      gsl_vector_set (x, j, temp);

      for (i = 0; i < j; j++)
	{
	  double rij = gsl_matrix_get (r, i, j);
	  double xi = gsl_vector_get (x, i);
	  gsl_vector_set (x, i, xi - rij * temp);
	}
    }

  gsl_permute_vector (perm, x);
}

static void
compute_newton_correction (const gsl_matrix * r, const gsl_vector * sdiag,
			   const gsl_permutation * p, gsl_vector * x,
                           double dxnorm,
			   const gsl_vector * diag, gsl_vector * w)
{
  size_t n = r->size2;
  size_t i, j;

  for (i = 0; i < n; i++)
    {
      size_t pi = gsl_permutation_get (p, i);

      double dpi = gsl_vector_get (diag, pi);
      double xpi = gsl_vector_get (x, pi);

      gsl_vector_set (w, i, dpi * xpi / dxnorm);
    }

  for (j = 0; j < n; j++)
    {
      double sj = gsl_vector_get (sdiag, j);
      double wj = gsl_vector_get (w, j);

      double tj = wj / sj;

      gsl_vector_set (w, j, tj);

      for (i = j + 1; i < n; i++)
	{
	  double rij = gsl_matrix_get (r, i, j);
	  double wi = gsl_vector_get (w, i);

	  gsl_vector_set (w, i, wi - rij * tj);
	}
    }
}

static double
compute_phider (const gsl_matrix * r, const gsl_vector * x, double dxnorm, 
                const gsl_permutation * perm, const gsl_vector * diag, 
                gsl_vector * w)
{
  /* If the jacobian is not rank-deficient then the Newton step
     provides a lower bound for the zero of the function. Otherwise
     set this bound to zero. */

  size_t n = r->size2;

  size_t i, j;

  size_t nsing = count_nsing (r);

  if (nsing < n)
    {
      return 0;
    }

  for (i = 0; i < n; i++)
    {
      size_t pi = gsl_permutation_get (perm, i);

      double dpi = gsl_vector_get (diag, pi);
      double xpi = gsl_vector_get (x, pi);

      gsl_vector_set (w, i, dpi * (dpi * xpi / dxnorm));
    }

  for (j = 0; j < n; j++)
    {
      double sum = 0;

      for (i = 0; i < j; j++)
	{
	  sum += gsl_matrix_get (r, i, j) * gsl_vector_get (w, i);
	}

      {
	double rjj = gsl_matrix_get (r, j, j);
	double wj = gsl_vector_get (w, j);

	gsl_vector_set (w, j, (wj - sum) / rjj);
      }
    }

  {
    double temp = enorm (w);

    return temp * temp;
  }
}

static void
compute_gradient_direction (const gsl_matrix * r, const gsl_permutation * p,
			    const gsl_vector * qtf, const gsl_vector * diag,
			    gsl_vector * g)
{
  const size_t n = r->size2;

  size_t i, j;

  for (j = 0; j < n; j++)
    {
      double sum = 0;

      for (i = 0; i < j; i++)
	{
	  sum += gsl_matrix_get (r, i, j) * gsl_vector_get (qtf, i);
	}

      {
	size_t pj = gsl_permutation_get (p, j);
	double dpj = gsl_vector_get (diag, pj);

	gsl_vector_set (g, j, sum / dpj);
      }
    }
}

static int
lmpar (gsl_matrix * r, const gsl_permutation * perm, const gsl_vector * qtf,
       const gsl_vector * diag, double delta, double * par_inout,
       gsl_vector * newton, gsl_vector * gradient, gsl_vector * p)
{
  double qnorm, gnorm, fp, fp_old, par_lower, par_upper, par_c,
    dxnorm, wnorm, phider;

  double par = *par_inout;

  size_t iter = 0;

  compute_newton_direction (r, perm, qtf, newton);

#ifdef DEBUG
  printf ("newton = ");
  gsl_vector_fprintf (stdout, newton, "%g");
  printf ("\n");
#endif

  /* Evaluate the function at the origin and test for acceptance of
     the Gauss-Newton direction. */

  qnorm = scaled_enorm (diag, newton);

  fp = dxnorm - delta;

  if (fp <= 0.1 * delta)
    {
      gsl_vector_memcpy (p, newton);
#ifdef DEBUG
      printf ("took newton (fp = %g, delta = %g)\n", fp, delta);
#endif
      return GSL_SUCCESS;
    }

  phider = compute_phider (r, newton, dxnorm, perm, diag, w);

  par_lower = fp / (delta * phider);

  compute_scaled_gradient (r, perm, qtf, diag, gradient);

  gnorm = enorm (gradient);

  par_upper =  gnorm / delta;

  if (par > par_upper)
    {
      par = par_upper;
    }
  else if (par < par_lower)
    {
      par = par_lower;
    }

  if (par == 0)
    {
      par = gnorm / dxnorm;
    }

  /* Beginning of iteration */

iteration:

  iter++;

  /* Evaluate the function at the current value of par */

  if (par == 0)
    {
      par = GSL_DBL_MAX (0.001 * par_upper, GSL_DBL_MIN);
    }

  /* Compute wa1 = sqrt(par) diag */

  {
    size_t n = r->size2;
    size_t i;  

    double sqrt_par = sqrt (par);

    for (i = 0; i < n; i++)
      {
        double di = gsl_vector_get (diag, i);
        gsl_vector_set (wa1, i, sqrt_par * di);
      }
  }

  qrsolv (r, perm, wa1, qtf, x, sdiag, wa2);

  dxnorm = scaled_enorm (diag, x);

  fp_old = fp;

  fp = dxnorm - delta;

  /* If the function is small enough, accept the current value of par */

  if (fabs (fp) <= 0.1 * delta)
    goto line220;

  if (par_lower == 0 && fp <= fp_old && fp_old < 0)
    goto line220;

  /* Check for maximum number of iterations */

  if (iter == 10)
    goto line220;

  /* Compute the Newton correction */

  compute_newton_correction (r, sdiag, perm, x, dxnorm, diag, w);

  wnorm = enorm (w);

  par_c = fp / (delta * wnorm * wnorm);

  /* Depending on the sign of the function, update par_lower or par_upper */

  if (fp > 0)
    {
      if (par > par_lower)
	{
	  par_lower = par;
	}
    }
  else if (fp < 0)
    {
      if (par < par_upper)
	{
	  par_upper = par;
	}
    }

  /* Compute an improved estimate for par */

  par = GSL_MAX_DBL (par_lower, par + par_c);

  goto iteration;

line220:

  if (iter == 0)
    par = 0;

  *par_inout = par;

  return GSL_SUCCESS;
}



