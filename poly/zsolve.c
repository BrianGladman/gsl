/* zsolve.c - finds the complex roots of  = 0 */

#include <config.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_poly.h>

#include <stdio.h>

#define INDEX(m,i,j,n) ((m)[(i)*(n) + (j)])

static void set_companion_matrix (const double *a, size_t n, double *m);
static void balance_companion_matrix (double *m, size_t n);

int
gsl_poly_complex_solve (const double *a, size_t n,
			gsl_poly_complex_workspace * w,
			gsl_complex_packed_ptr z)
{
  double *m;

  m = w->matrix;

  set_companion_matrix (a, n - 1, m);

  dump (w);

  balance_companion_matrix (m, n - 1);

  dump (w);
}

static void
set_companion_matrix (const double *a, size_t nc, double *m)
{
  size_t i, j;

  for (i = 0; i < nc; i++)
    for (j = 0; j < nc; j++)
      INDEX (m, i, j, nc) = 0.0;

  for (i = 1; i < nc; i++)
    INDEX (m, i, i - 1, nc) = 1.0;

  for (i = 0; i < nc; i++)
    INDEX (m, i, nc - 1, nc) = -a[i] / a[nc];
}

#include <stdio.h>

#define RADIX 2
#define RADIX2 (RADIX*RADIX)

static void
balance_companion_matrix (double *m, size_t nc)
{
  int not_converged = 1;

  double row_norm = 0;
  double column_norm = 0;

  while (not_converged)
    {
      size_t i, j;
      double g, f, s;

      not_converged = 0;

      for (i = 0; i < nc; i++)
	{
	  size_t ii, jj;

	  printf ("i = %d\n", i);

	  for (ii = 0; ii < nc; ii++)
	    {
	      for (jj = 0; jj < nc; jj++)
		printf (" %10g", INDEX (m, ii, jj, nc));

	      printf ("\n");
	    }
	  printf ("\n");


	  /* column norm, excluding the diagonal */

	  if (i != nc - 1)
	    {
	      column_norm = fabs (INDEX (m, i + 1, i, nc));
	    }
	  else
	    {
	      column_norm = 0;

	      for (j = 0; j < nc - 1; j++)
		{
		  column_norm += fabs (INDEX (m, j, nc - 1, nc));
		}
	    }

	  /* row norm, excluding the diagonal */

	  if (i == 0)
	    {
	      row_norm = fabs (INDEX (m, 0, nc - 1, nc));
	    }
	  else if (i == nc - 1)
	    {
	      row_norm = fabs (INDEX (m, i, i - 1, nc));
	    }
	  else
	    {
	      row_norm = fabs (INDEX (m, i, i - 1, nc)) + fabs (INDEX (m, i, nc - 1, nc));
	    }

	  if (column_norm == 0 || row_norm == 0)
	    {
	      continue;
	    }

	  g = row_norm / RADIX;
	  f = 1;
	  s = column_norm + row_norm;

	  while (column_norm < g)
	    {
	      f *= RADIX;
	      column_norm *= RADIX2;
	    }

	  g = row_norm * RADIX;

	  while (column_norm > g)
	    {
	      f /= RADIX;
	      column_norm /= RADIX2;
	    }

	  if ((row_norm + column_norm) < 0.95 * s * f)
	    {
	      not_converged = 1;

	      g = 1 / f;

	      if (i == 0)
		{
		  INDEX (m, 0, nc - 1, nc) *= g;
		}
	      else
		{
		  INDEX (m, i, i - 1, nc) *= g;
		  INDEX (m, i, nc - 1, nc) *= g;
		}

	      if (i == nc - 1)
		{
		  for (j = 0; j < nc; j++)
		    {
		      INDEX (m, j, i, nc) *= f;
		    }
		}
	      else
		{
		  INDEX (m, i + 1, i, nc) *= f;
		}
	    }
	}
    }
}

static double
norm_companion (const double *m, size_t nc)
{
  size_t i;
  double A, s2 = 0;

  for (i = 1; i < nc - 1; i++)
    {
      const double t = INDEX (m, i + 1, i, nc);
      s2 += t * t;
    }

  for (i = 0; i < nc; i++)
    {
      const double t = INDEX (m, i, nc, nc);
      s2 += t * t;
    }

  A = sqrt (s2);

  return A;
}


static void
qr_companion (double *h, size_t nc, double *wr, double *wi)
{

  n = nc;
  t = 0.0;

nextw:

  if (n == 0)
    goto 9;

  its = 0;

  na = n - 1;

label2:

  for (k = n - 1; k != 0; k--)
    {
      double a1 = fabs (INDEX (h, k, k - 1, nc));
      double a2 = fabs (INDEX (h, k - 1, k - 1, nc));
      double a3 = fabs (INDEX (h, k, k, nc));

      if (a1 <= GSL_DBL_EPSILON * (a2 + a3))
	break;
    }

  x = INDEX (h, n, n, nc);

  if (k == n)
    {
      wr[n] = x + t;		/* one real root */
      wi[n] = 0.0;
      n = na;
      continue;
    }

  y = INDEX (h, na, na, nc);

  w = INDEX (h, n, na, nc) * INDEX (h, na, n, nc);

  if (k == na)
    {
      p = (y - x) / 2;
      q = p * p + w;
      y = sqrt (fabs (q));

      x += t;

      if (q > 0)		/* two real roots */
	{
	  if (p < 0)
	    y = -y;
	  y += p;

	  wr[n] = x - w / y;
	  wi[n] = 0;
	  wr[na] = x + y;
	  wi[na] = 0;
	}
      else
	{
	  wr[n] = x + p;
	  wr[na] = x + p;
	  wi[n] = -y;
	  wi[na] = y;
	}
      n = n - 2;

      continue;
    }

  /* No more roots found yet, do another iteration */

  if (its == 10 || its == 20)
    {
      /* use an exceptional shift */

      t += x;

      for (i = 0; i < n; i++)
	{
	  INDEX (h, i, i, nc) -= x;
	}

      s = fabs (INDEX (h, n, na, nc)) + fabs (INDEX (h, na, n - 2, nc));
      y = 0.75 * s;
      x = y;
      w = -0.4375 * s * s;
    }

  its++;

  for (m = n - 2; m >= k; m--)
    {
      z = INDEX (h, m, m, nc);
      r = x - z;
      s = y - z;
      p = INDEX (h, m, m + 1, nc) + (r * s - w) / INDEX (h, m + 1, m, nc);
      q = INDEX (h, m + 1, m + 1, nc) - z - r - s;
      r = INDEX (h, m + 2, m + 1, nc);
      s = fabs (p) + fabs (q) + fabs (r);
      p /= s;
      q /= s;
      r /= s;

      if (m == k)
	break;

      a1 = fabs(INDEX(h, m, m - 1, nc)) ;
      a2 = fabs(INDEX(h,m-1,m-1,nc)) + fabs(INDEX(h,m+1,m+1,nc));

      if (a1 * (fabs(q) + fabs(r)) <= GSL_DBL_EPSILON * fabs(p) * a2)
        break ;
    }

  for (i = m + 2; i < n ; i++)
    {
      INDEX(h,i,i-2,nc) = 0;
    }

  for (i = m + 3; i < n ; i++)
    {
      INDEX(h,i,i-3,nc) = 0;
    }

  /* double QR step */

