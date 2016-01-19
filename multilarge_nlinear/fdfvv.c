/* multilarge_nlinear/fdfvv.c
 * 
 * Copyright (C) 2015 Patrick Alken
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

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/*
fdfvv()
  Compute approximate second directional derivative using
finite differences, specifically J(x)^T f_vv(x)

See Eq. 19 of:

M. K. Transtrum, J. P. Sethna, Improvements to the Levenberg
Marquardt algorithm for nonlinear least-squares minimization,
arXiv:1201.5885, 2012.

Inputs: h     - step size for finite difference
        x     - parameter vector, size p
        v     - geodesic velocity, size p
        g     - gradient J(x)^T f(x), size p
        JTJ   - J(x)^T J(x), p-by-p
        fdf   - fdf struct
        JTfvv - (output) J(x)^T f_vv(x)
        workn - workspace, size n
        workp - workspace, size p

Return: success or error
*/

static int
fdfvv(const double h, const gsl_vector *x, const gsl_vector *v,
      const gsl_vector *g, const gsl_matrix *JTJ,
      gsl_multilarge_nlinear_fdf *fdf,
      gsl_vector *JTfvv, gsl_vector *workn, gsl_vector *workp)
{
  int status;
  const size_t p = fdf->p;
  const double hinv = 1.0 / h;
  gsl_vector_view JTJv = gsl_vector_subvector(workn, 0, p);
  size_t i;

  /* compute workp = x + h*v */
  for (i = 0; i < p; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double vi = gsl_vector_get(v, i);

      gsl_vector_set(workp, i, xi + h * vi);
    }

  /* compute workn = f(x + h*v) */
  status = gsl_multilarge_nlinear_eval_f (fdf, workp, workn);
  if (status)
    return status;

  /* compute workp = J(x)^T f(x + h*v) */
  status = gsl_multilarge_nlinear_eval_df(fdf, x, workn, workp, NULL);
  if (status)
    return status;

  /* compute JTJv = J^T J v */
  gsl_blas_dsymv(CblasLower, 1.0, JTJ, v, 0.0, &JTJv.vector);

  for (i = 0; i < p; ++i)
    {
      double JTfpi = gsl_vector_get(workp, i);        /* J^T f(x + h*v) */
      double gi = gsl_vector_get(g, i);               /* J^T f(x) */
      double JTJvi = gsl_vector_get(&JTJv.vector, i); /* J^T J v */
      double u;

      u = (2.0 * hinv) * ((JTfpi - gi) * hinv - JTJvi);
      gsl_vector_set(JTfvv, i, u);
    }

  return GSL_SUCCESS;
}

/*
gsl_multilarge_nlinear_fdJTfvv()
  Compute approximate second directional derivative
vector J(x)^T f_vv(x) using finite differences

Inputs: h     - step size for finite difference
        x     - parameter vector, size p
        v     - geodesic velocity, size p
        g     - gradient J(x)^T f(x), size p
        JTJ   - J(x)^T J(x), p-by-p
        fdf   - fdf
        JTfvv - (output) approximate J(x)^T f_vv(x), size p
        workn - workspace, size n
        workp - workspace, size p

Return: success or error
*/

int
gsl_multilarge_nlinear_fdJTfvv(const double h, const gsl_vector *x, const gsl_vector *v,
                               const gsl_vector *g, const gsl_matrix *JTJ,
                               gsl_multilarge_nlinear_fdf *fdf,
                               gsl_vector *JTfvv, gsl_vector *workn, gsl_vector *workp)
{
  return fdfvv(h, x, v, g, JTJ, fdf, JTfvv, workn, workp);
}
