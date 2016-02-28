/* multifit_nlinear/common.c
 * 
 * Copyright (C) 2014, 2015, 2016 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

static double scaled_norm(const gsl_vector *D, const gsl_vector *a);
static void scaled_addition (const double alpha, const gsl_vector * x,
                             const double beta, const gsl_vector * y,
                             gsl_vector * z);

/* compute || diag(D) a || */
static double
scaled_norm(const gsl_vector *D, const gsl_vector *a)
{
  const size_t n = a->size;
  double e2 = 0.0;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double Di = gsl_vector_get(D, i);
      double ai = gsl_vector_get(a, i);
      double u = Di * ai;

      e2 += u * u;
    }

  return sqrt (e2);
}

/* compute z = alpha*x + beta*y */
static void
scaled_addition (const double alpha, const gsl_vector * x,
                 const double beta, const gsl_vector * y, gsl_vector * z)
{
  const size_t N = z->size;
  size_t i;

  for (i = 0; i < N; i++)
    {
      double xi = gsl_vector_get (x, i);
      double yi = gsl_vector_get (y, i);
      gsl_vector_set (z, i, alpha * xi + beta * yi);
    }
}
