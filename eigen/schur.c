/* eigen/schur.c
 * 
 * Copyright (C) 2006, 2007 Patrick Alken
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * This module contains some routines related to manipulating the
 * Schur form of a matrix which are needed by the eigenvalue solvers
 *
 * This file contains routines based on original code from LAPACK
 * which is distributed under the modified BSD license. The LAPACK
 * routine used is DLANV2.
 */

static inline void schur_standard_form(gsl_matrix *A, double *cs, double *sn);

/*******************************************************
 *            INTERNAL ROUTINES                        *
 *******************************************************/

/*
schur_standard_form()
  Compute the Schur factorization of a real 2-by-2 matrix in
standard form:

[ A B ] = [ CS -SN ] [ T11 T12 ] [ CS SN ]
[ C D ]   [ SN  CS ] [ T21 T22 ] [-SN CS ]

where either:
1) T21 = 0 so that T11 and T22 are real eigenvalues of the matrix, or
2) T11 = T22 and T21*T12 < 0, so that T11 +/- sqrt(|T21*T12|) are
   complex conjugate eigenvalues

Inputs: A     - 2-by-2 matrix
        cs    - where to store cosine parameter of rotation matrix
        sn    - where to store sine parameter of rotation matrix

Notes: 1) based on LAPACK routine DLANV2
       2) On output, A is modified to contain the matrix in standard form
*/

static inline void
schur_standard_form(gsl_matrix *A, double *cs, double *sn)
{
  double a, b, c, d; /* input matrix values */
  double tmp;
  double p, z;
  double bcmax, bcmis, scale;
  double tau, sigma;
  double cs1, sn1;
  double aa, bb, cc, dd;
  double sab, sac;

  a = gsl_matrix_get(A, 0, 0);
  b = gsl_matrix_get(A, 0, 1);
  c = gsl_matrix_get(A, 1, 0);
  d = gsl_matrix_get(A, 1, 1);

  if (c == 0.0)
    {
      /*
       * matrix is already upper triangular - set rotation matrix
       * to the identity
       */
      *cs = 1.0;
      *sn = 0.0;
    }
  else if (b == 0.0)
    {
      /* swap rows and columns to make it upper triangular */

      *cs = 0.0;
      *sn = 1.0;

      tmp = d;
      d = a;
      a = tmp;
      b = -c;
      c = 0.0;
    }
  else if (((a - d) == 0.0) && (GSL_SIGN(b) != GSL_SIGN(c)))
    {
      /* the matrix has complex eigenvalues with a == d */
      *cs = 1.0;
      *sn = 0.0;
    }
  else
    {
      tmp = a - d;
      p = 0.5 * tmp;
      bcmax = GSL_MAX(fabs(b), fabs(c));
      bcmis = GSL_MIN(fabs(b), fabs(c)) * GSL_SIGN(b) * GSL_SIGN(c);
      scale = GSL_MAX(fabs(p), bcmax);
      z = (p / scale) * p + (bcmax / scale) * bcmis;

      if (z >= 4.0 * GSL_DBL_EPSILON)
        {
          /* real eigenvalues, compute a and d */

          z = p + GSL_SIGN(p) * fabs(sqrt(scale) * sqrt(z));
          a = d + z;
          d -= (bcmax / z) * bcmis;

          /* compute b and the rotation matrix */

          tau = gsl_hypot(c, z);
          *cs = z / tau;
          *sn = c / tau;
          b -= c;
          c = 0.0;
        }
      else
        {
          /*
           * complex eigenvalues, or real (almost) equal eigenvalues -
           * make diagonal elements equal
           */

          sigma = b + c;
          tau = gsl_hypot(sigma, tmp);
          *cs = sqrt(0.5 * (1.0 + fabs(sigma) / tau));
          *sn = -(p / (tau * (*cs))) * GSL_SIGN(sigma);

          /*
           * Compute [ AA BB ] = [ A B ] [ CS -SN ]
           *         [ CC DD ]   [ C D ] [ SN  CS ]
           */
          aa = a * (*cs) + b * (*sn);
          bb = -a * (*sn) + b * (*cs);
          cc = c * (*cs) + d * (*sn);
          dd = -c * (*sn) + d * (*cs);

          /*
           * Compute [ A B ] = [ CS SN ] [ AA BB ]
           *         [ C D ]   [-SN CS ] [ CC DD ]
           */
          a = aa * (*cs) + cc * (*sn);
          b = bb * (*cs) + dd * (*sn);
          c = -aa * (*sn) + cc * (*cs);
          d = -bb * (*sn) + dd * (*cs);

          tmp = 0.5 * (a + d);
          a = d = tmp;

          if (c != 0.0)
            {
              if (b != 0.0)
                {
                  if (GSL_SIGN(b) == GSL_SIGN(c))
                    {
                      /*
                       * real eigenvalues: reduce to upper triangular
                       * form
                       */
                      sab = sqrt(fabs(b));
                      sac = sqrt(fabs(c));
                      p = GSL_SIGN(c) * fabs(sab * sac);
                      tau = 1.0 / sqrt(fabs(b + c));
                      a = tmp + p;
                      d = tmp - p;
                      b -= c;
                      c = 0.0;

                      cs1 = sab * tau;
                      sn1 = sac * tau;
                      tmp = (*cs) * cs1 - (*sn) * sn1;
                      *sn = (*cs) * sn1 + (*sn) * cs1;
                      *cs = tmp;
                    }
                }
              else
                {
                  b = -c;
                  c = 0.0;
                  tmp = *cs;
                  *cs = -(*sn);
                  *sn = tmp;
                }
            }
        }
    }

  /* set new matrix elements */

  gsl_matrix_set(A, 0, 0, a);
  gsl_matrix_set(A, 0, 1, b);
  gsl_matrix_set(A, 1, 0, c);
  gsl_matrix_set(A, 1, 1, d);
} /* schur_standard_form() */
