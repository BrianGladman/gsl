/* poly/oper.c
 * 
 * Copyright (C) 2002 Gert Van den Eynde
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>

extern int gsl_check_range;	/* defined in vector/vector.c */

double
gsl_poly_eval2 (gsl_poly * p, double x)
{
  return gsl_poly_eval (p->c, p->size, x);
}

double
gsl_poly_get (const gsl_poly * p, const size_t i)
{
  if (gsl_check_range)
    {
      if (i >= p->size)
	{
	  GSL_ERROR_VAL ("index out of range", GSL_EINVAL, 0.0);
	}
    }

  return p->c[i];
}

void
gsl_poly_set (gsl_poly * p, const size_t i, double a)
{
  if (gsl_check_range)
    {
      if (i >= p->size)
	{
	  GSL_ERROR_VOID ("index out of range", GSL_EINVAL);
	}
    }

  p->c[i] = a;
}

int
gsl_poly_set_zero (gsl_poly * p)
{
  size_t i;
  size_t size = p->size;

  for (i = 0; i < size; i++)
    {
      p->c[i] = 0.0;
    }

  return GSL_SUCCESS;
}

int
gsl_poly_memcpy (gsl_poly * dest, const gsl_poly * src)
{
  size_t src_size = src->size;
  size_t dest_size = dest->size;
  size_t j;

  if (src_size != dest_size)
    {
      GSL_ERROR ("polynomial sizes are not equal", GSL_EBADLEN);
    }

  for (j = 0; j < src_size; j++)
    {
      dest->c[j] = src->c[j];
    }

  return GSL_SUCCESS;
}

int
gsl_poly_chop (gsl_poly * p, size_t k)
{
  size_t i, size = p->size;

  for (i = k; i < size; i++)
    {
      p->c[i] = 0.0;
    }

  return GSL_SUCCESS;
}

size_t gsl_poly_find_size (const gsl_poly * p, double tol)
{
  size_t i;
  size_t size = p->size;

  i = size;

  while (i > 0 && fabs (p->c[i-1]) < tol)
    {
      i--;
    }

  return i;
}

int
gsl_poly_set_from_array (gsl_poly * p, const double a[], size_t size)
{
  size_t i;

  if (p->size != size)
    {
      GSL_ERROR ("size of polynomial must match array", GSL_EBADLEN);
    }

  for (i = 0; i < size; i++)
    {
      p->c[i] = a[i];
    }

  return GSL_SUCCESS;
}

int
gsl_poly_add (gsl_poly * p1, const gsl_poly * p2)
{
  size_t i;
  size_t size = p2->size;

  if (p1->size != p2->size)
    {
      GSL_ERROR ("size of polynomials must be equal", GSL_EBADLEN);
    }

  for (i = 0; i < size; i++)
    {
      p1->c[i] += p2->c[i];
    }

  return GSL_SUCCESS;
}

int
gsl_poly_sub (gsl_poly * p1, const gsl_poly * p2)
{
  size_t i;
  size_t size = p2->size;

  if (p1->size != p2->size)
    {
      GSL_ERROR ("size of polynomials must be equal", GSL_EBADLEN);
    }

  for (i = 0; i < size; i++)
    {
      p1->c[i] -= p2->c[i];
    }

  return GSL_SUCCESS;
}

int
gsl_poly_mul (gsl_poly * q, const gsl_poly * p1, const gsl_poly * p2)
{
  int i, j;

  int size_p1 = p1->size;
  int size_p2 = p2->size;
  int size_q = q->size;

  if (size_q + 1 != size_p1 + size_p2)
    {
      GSL_ERROR ("incorrect size of destination polynomial", GSL_EBADLEN);
    }

  for (i = 0; i < size_q; i++)
    {
      q->c[i] = 0.0;
    }

  for (i = 0; i < size_q; i++)
    {
      double sum = 0.0;

      size_t j_min = (i >= size_p2) ? (i - size_p2 + 1) : 0;
      size_t j_max = (i >= size_p1) ? size_p1 : (i + 1);

      for (j = j_min; j < j_max; j++)
	{
	  sum += p1->c[j] * p2->c[i - j];
	}

      q->c[i] = sum;
    }

  return GSL_SUCCESS;
}

int
gsl_poly_div (gsl_poly * q, gsl_poly * r, const gsl_poly * u,
	      const gsl_poly * v)
{
  size_t size_q = q->size;
  size_t size_r = r->size;
  size_t size_u = u->size;
  size_t size_v = v->size;

  size_t i, j, k;

  if (size_v > size_u)
    {
      GSL_ERROR
	("degree of denominator should not exceed degree of numerator",
	 GSL_EINVAL);
    }

  if (size_q != size_u)
    {
      GSL_ERROR ("size of quotient must equal size of numerator",
		 GSL_EBADLEN);
    }

  if (size_r != size_u)
    {
      GSL_ERROR ("size of remainder must equal size of numerator",
		 GSL_EBADLEN);
    }

  /* Trim leading zeros from denominator */

  while (size_v > 1 && (v->c[size_v - 1] == 0.0))
    {
      size_v--;
    }

  for (i = 0; i < size_u; i++)
    {
      q->c[i] = 0.0;
      r->c[i] = u->c[i];
    }

  for (k = size_u; k >= size_v && k--;)
    {
      double qk = r->c[k] / v->c[size_v - 1];

      q->c[k - (size_v - 1)] = qk;

      /* Shift the remainder */

      r->c[k] = 0.0;

      for (j = 0; j < size_v - 1; j++)
	{
	  r->c[j + (k - (size_v - 1))] -= qk * v->c[j];
	}
    }

  return GSL_SUCCESS;
}

int
gsl_poly_scale (gsl_poly * p, double a)
{
  size_t size = p->size;
  size_t i;

  for (i = 0; i < size; i++)
    {
      p->c[i] *= a;
    }

  return GSL_SUCCESS;
}

int
gsl_poly_diff (gsl_poly * dp, const gsl_poly * p)
{
  size_t size = dp->size;
  size_t i;

  if (size + 1 != p->size)
    {
      GSL_ERROR ("derivative must be one element smaller than polynomial",
		 GSL_EBADLEN);
    }

  for (i = 0; i < size; i++)
    {
      dp->c[i] = (i + 1.0) * p->c[i + 1];
    }

  return GSL_SUCCESS;
}

void
gsl_poly_dump (gsl_poly * p)
{
  int i;
  printf ("%d\n", p->size);
  for (i = 0; i < p->size; i++)
    {
      printf ("%e\n", p->c[i]);
    }
  printf ("\n");
}
