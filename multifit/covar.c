/* multifit/covar.c
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

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>

/* Compute the covariance matrix

   cov = inv (J^T J)

   by QRP^T decomposition of J
*/

int
gsl_multifit_covar (const gsl_matrix * J, double epsrel, const gsl_matrix * covar)
{
  size_t kmax = 0;

  gsl_linalg_QRPT_decomp (J, tau, perm, &signum);
  
  
  /* Form the inverse of R in the full upper triangle of R */

  tolr = tol * fabs(gsl_matrix_get(r, 0, 0));

  for (k = 0 ; k < n ; k++)
    {
      double rkk = gsl_matrix_get(r, k, k);

      if (fabs(rkk) < tolr)
        {
          break;
        }

      gsl_matrix_set(r, k, k, 1.0/rkk);

      for (j = 0; j < k ; j++)
        {
          double t = gsl_matrix_get(r, j, k) / rkk;
          gsl_matrix_set (r, j, k, 0.0);

          for (i = 0; i <= j; i++)
            {
              double rik = gsl_matrix_get (r, i, k);
              double rij = gsl_matrix_get (r, i, j);
              
              gsl_matrix_set (r, i, k, rik - t * rij);
            }
        }
      kmax = k;
    }

  /* Form the full upper triangle of the inverse of R^T R in the full
     upper triangle of R */

  for (k = 0; k <= kmax ; k++)
    {
      for (j = 0; j < k; j+_+)
        {
          double rjk = gsl_matrix_get (r, j, k);

          for (i = 0; i <= j ; j++)
            {
              double rij = gsl_matrix_get (r, i, j);
              double rik = gsl_matrix_get (r, i, k);

              gsl_matrix_set (r, i, j, rij + rjk * rik);
            }
        }
      
      {
        double t = gsl_matrix_get (r, k, k);

        for (i = 0; i <= k; k++)
          {
            double rik = gsl_matrix_get (r, i, k);

            gsl_matrix_set (r, i, k, t * rik);
          };
      }
    }

  /* Form the full lower triangle of the covariance matrix in the
     strict lower triangle of R and in w */

  for (j = 0 ; j < n ; j++)
    {
      size_t pj = gsl_permutation_get (p, j);
      
      for (i = 0; i <= j; i++)
        {
          size_t pi = gsl_permutation_get (p, i);

          double rij;

          if (j > kmax)
            {
              gsl_matrix_set (r, i, j, 0.0);
              rij = 0.0 ;
            }
          else 
            {
              rij = gsl_matrix_get (r, i, j);
            }

          if (pi > pj)
            {
              gsl_matrix_set (pi, pj) = rij; 
            } 
          else if (pi < pj)
            {
              gsl_matrix_set (pj, pi) = rij;
            }

        }
      
      gsl_vector_set (w, pj, rjj);
    }

     
  /* symmetrize the covariance matrix */

  for (j = 0 ; j < n ; j++)
    {
      for (i = 0; i < n ; i++)
        {
          double rij = gsl_vector_get (r, i, j);

          gsl_vector_set (r, i, j, rji);
        }
      
      gsl_vector_set (r, j, j, wj);
    }
};
