/* legendre_P.c
 * 
 * Copyright (C) 2009-2013 Patrick Alken
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

/*
 * The routines in this module compute associated Legendre functions
 * (ALFs) up to order and degree 2700, using the method described
 * in
 *
 * [1] S. A. Holmes and W. E. Featherstone, A unified approach
 *     to the Clenshaw summation and the recursive computation of very
 *     high degree and order normalised associated Legendre functions,
 *     Journal of Geodesy, 76, pg. 279-299, 2002.
 *
 * Further information on ALFs can be found in
 *
 * [2] Abramowitz and Stegun, Handbook of Mathematical Functions,
 *     Chapter 8, 1972.
 */

static void legendre_sqrts(const size_t lmax, double *array);

/* number of P_{lm} functions for a given lmax */
size_t
gsl_sf_legendre_nlm(const size_t lmax)
{
  return ((lmax + 1) * (lmax + 2) / 2);
}

/*
gsl_sf_legendre_array_n()
  This routine returns the minimum result_array[] size needed
for a given lmax
*/

size_t
gsl_sf_legendre_array_n(const size_t lmax)
{
  size_t nlm = gsl_sf_legendre_nlm(lmax);
  size_t nsqrt = 2 * lmax + 2; /* extra room to precompute sqrt factors */

  return (nlm + nsqrt);
} /* gsl_sf_legendre_array_n() */

/*
gsl_sf_legendre_array_index()
This routine computes the index into a result_array[] corresponding
to a given (l,m)
*/

size_t
gsl_sf_legendre_array_index(const size_t l, const size_t m)
{
  return (l * (l + 1) / 2 + m);
} /* alf_index() */

/*
gsl_sf_legendre_array_schmidt()
  Compute Schmidt semi-normalized ALFs

Inputs: lmax         - maximum degree
        x            - input argument in [-1,1]
        result_array - (output) where to store S_{lm}(x)
*/

int
gsl_sf_legendre_array_schmidt(const size_t lmax, const double x,
                              double result_array[])
{
  int s = gsl_sf_legendre_array_schmidt_e(lmax, x, 1.0, result_array);

  return s;
} /* gsl_sf_legendre_array_schmidt() */

/*
gsl_sf_legendre_deriv_array_schmidt()
  Compute Schmidt semi-normalized ALFs and their first derivatives
S_{nm}(x) and S'_{nm}(x)

Inputs: lmax               - maximum degree
        x                  - input argument
        result_array       - (output) where to store S_{lm}(x)
        result_deriv_array - (output) where to store d/dx S_{lm}(x)
*/

int
gsl_sf_legendre_deriv_array_schmidt(const size_t lmax, const double x,
                                    double result_array[],
                                    double result_deriv_array[])
{
  int s;
  const double u = sqrt((1.0 - x) * (1.0 + x));
  const double uinv = 1.0 / u;
  const size_t n = gsl_sf_legendre_array_n(lmax);
  size_t i;

  /* compute d/dtheta S_{lm}(x) */
  s = gsl_sf_legendre_deriv_alt_array_schmidt_e(lmax, x, 1.0, result_array,
                                                result_deriv_array);

  /* scale derivative array by u to recover P'(x) from d/dtheta P(x) */
  for (i = 0; i < n; ++i)
    result_deriv_array[i] *= -uinv;

  return s;
} /* gsl_sf_legendre_deriv_array_schmidt() */

/*
gsl_sf_legendre_deriv2_array_schmidt()
  Compute Schmidt semi-normalized ALFs and their first derivatives
S_{nm}(x) and S'_{nm}(x)

Inputs: lmax                - maximum degree
        x                   - input argument
        result_array        - (output) where to store S_{lm}(x)
        result_deriv_array  - (output) where to store d/dx S_{lm}(x)
        result_deriv2_array - (output) where to store d^2/dx^2 S_{lm}(x)
*/

int
gsl_sf_legendre_deriv2_array_schmidt(const size_t lmax, const double x,
                                     double result_array[],
                                     double result_deriv_array[],
                                     double result_deriv2_array[])
{
  int s;
  const double u = sqrt((1.0 - x) * (1.0 + x));
  const double uinv = 1.0 / u;
  const double uinv2 = uinv * uinv;
  const size_t n = gsl_sf_legendre_array_n(lmax);
  size_t i;

  /* compute d/dtheta S_{lm}(x) and d^2/dtheta^2 S_{lm}(x) */
  s = gsl_sf_legendre_deriv2_alt_array_schmidt_e(lmax, x, 1.0, result_array,
                                                 result_deriv_array,
                                                 result_deriv2_array);
  if (s)
    return s;

  /* scale derivative arrays to recover P'(x) and P''(x) */
  for (i = 0; i < n; ++i)
    {
      double dp = result_deriv_array[i];
      double d2p = result_deriv2_array[i];

      result_deriv_array[i] *= -uinv;
      result_deriv2_array[i] = (d2p - x * uinv * dp) * uinv2;
    }

  return s;
} /* gsl_sf_legendre_deriv2_array_schmidt() */

int
gsl_sf_legendre_array_spharm(const size_t lmax, const double x,
                             double result_array[])
{
  int s = gsl_sf_legendre_array_spharm_e(lmax, x, -1.0, result_array);

  return s;
}

/*
gsl_sf_legendre_array_spharm_e()
  Compute spherical harmonic normalized ALFs

Inputs: lmax         - maximum degree
        x            - input argument in [-1,1]
        csphase      - include Condon-Shortley phase?
        result_array - (output) where to store Y_{lm}(x)
*/

int
gsl_sf_legendre_array_spharm_e(const size_t lmax, const double x,
                               const double csphase,
                               double result_array[])
{
  int s;
  const double fac1 = 1.0 / sqrt(4.0 * M_PI);
  const double fac2 = 1.0 / sqrt(8.0 * M_PI);
  size_t l, m;
  size_t twoellp1; /* 2l + 1 */
  size_t nlm = gsl_sf_legendre_nlm(lmax);
  double *sqrts = &(result_array[nlm]);
  
  s = gsl_sf_legendre_array_schmidt_e(lmax, x, csphase, result_array);
  if (s)
    return s;

  twoellp1 = 1;
  for (l = 0; l <= lmax; ++l)
    {
      /*
       * handle m = 0 case separately since S_{lm} have a
       * different normalization for m = 0 while Y_{lm} do not
       */
      result_array[gsl_sf_legendre_array_index(l, 0)] *= sqrts[twoellp1] * fac1;

      for (m = 1; m <= l; ++m)
        result_array[gsl_sf_legendre_array_index(l, m)] *= sqrts[twoellp1] * fac2;

      twoellp1 += 2;
    }

  return s;
} /* gsl_sf_legendre_array_spharm_e() */

/*********************************************************
 *                 INTERNAL ROUTINES                     *
 *********************************************************/

#if 0
/*
gsl_Plm_array_none()
  This routine computes unnormalized associated Legendre polynomials
*/

static int
gsl_Plm_array_none(const size_t lmax, const double x,
                   double result_array[])
{
  const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
  size_t l, m;
  size_t k, kstart;
  double plm, pmm;
  double pm1,    /* P(l-1,m) */
         pm2;    /* P(l-2,m) */
  double twomm1; /* 2*m - 1 */

  /* initial values P(0,0) and P(1,0) */

  pm2 = 1.0; /* P(0,0) */
  pm1 = x;   /* P(1,0) */

  result_array[0] = pm2;

  if (lmax == 0)
    return 0;

  result_array[1] = pm1;

  /* Compute P(l,0) */

  k = 1;
  for (l = 2; l <= lmax; ++l)
    {
      k += l;
      plm = ((2*l - 1) * x * pm1 - (l - 1) * pm2) / (double) l;
      result_array[k] = plm;
      pm2 = pm1;
      pm1 = plm;
    }

  /* Compute P(m,m), P(m+1,m) and P(l,m) */

  pmm = 1.0;
  twomm1 = -1.0; /* 2 * m - 1 */

  kstart = 0;

  for (m = 1; m <= lmax - 1; ++m)
    {
      /*
       * compute
       *
       * P(m,m) = u * (2m - 1) P(m-1,m-1)
       * and
       * dP(m,m)/dx = -m * x * P(m,m) / u^2
       */
      kstart += m + 1;
      twomm1 += 2.0;
      pmm *= w->csphase * u * twomm1;
      result_array[kstart] = pmm;
      pm2 = pmm;

      /*
       * compute
       *
       * P(m+1,m) = (2 * m + 1) * x * P(m,m)
       * and
       * dP(m+1,m)/dx = [(2*m + 1) * P(m,m) - (m+1) * x * P(m+1,m)]/u^2
       */
      k = kstart + m + 1;
      pm1 = x * pmm * (2*m + 1);
      result_array[k] = pm1;

      /* compute P(l,m) */
      for (l = m + 2; l <= lmax; ++l)
        {
          k += l;
          plm = ((2*l - 1) * x * pm1 - (l + m - 1) * pm2) /
                (double) (l - m);
          result_array[k] = plm;
          pm2 = pm1;
          pm1 = plm;
        }
    } /* for (m = 1; m <= lmax - 1; ++m) */

  /* compute P(lmax,lmax) */

  kstart += m + 1;
  twomm1 += 2.0;
  pmm *= w->csphase * u * twomm1;
  result_array[kstart] = pmm;

  return 0;
} /* gsl_Plm_array_none() */

/*
gsl_Plm_deriv_array_none()
  This routine computes unnormalized associated Legendre polynomials
and their first derivatives.
*/

static int
gsl_Plm_deriv_array_none(const size_t lmax, const double x,
                         double result_array[],
                         double result_deriv_array[], alf_workspace *w)
{
  const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
  const double u_sq = u * u;
  size_t l, m;
  size_t k, kstart;
  double plm, pmm;
  double pm1,    /* P(l-1,m) */
         pm2;    /* P(l-2,m) */
  double twomm1; /* 2*m - 1 */

  /* initial values P(0,0) and P(1,0) */

  pm2 = 1.0; /* P(0,0) */
  pm1 = x;   /* P(1,0) */

  result_array[0] = pm2;
  result_deriv_array[0] = 0.0;

  if (lmax == 0)
    return 0;

  result_array[1] = pm1;
  result_deriv_array[1] = 1.0;

  /* Compute P(l,0) */

  k = 1;
  for (l = 2; l <= lmax; ++l)
    {
      k += l;
      plm = ((2*l - 1) * x * pm1 - (l - 1) * pm2) / (double) l;
      result_array[k] = plm;
      result_deriv_array[k] = l * (pm1 - x * plm) / u_sq;
      pm2 = pm1;
      pm1 = plm;
    }

  /* Compute P(m,m), P(m+1,m) and P(l,m) */

  pmm = 1.0;
  twomm1 = -1.0; /* 2 * m - 1 */

  kstart = 0;

  for (m = 1; m <= lmax - 1; ++m)
    {
      /*
       * compute
       *
       * P(m,m) = u * (2m - 1) P(m-1,m-1)
       * and
       * dP(m,m)/dx = -m * x * P(m,m) / u^2
       */
      kstart += m + 1;
      twomm1 += 2.0;
      pmm *= w->csphase * u * twomm1;
      result_array[kstart] = pmm;
      result_deriv_array[kstart] = - m * x * pmm / u_sq;
      pm2 = pmm;

      /*
       * compute
       *
       * P(m+1,m) = (2 * m + 1) * x * P(m,m)
       * and
       * dP(m+1,m)/dx = [(2*m + 1) * P(m,m) - (m+1) * x * P(m+1,m)]/u^2
       */
      k = kstart + m + 1;
      pm1 = x * pmm * (2*m + 1);
      result_array[k] = pm1;
      result_deriv_array[k] = ((2*m + 1) * pmm - (m + 1) * x * pm1) / u_sq;

      /* compute P(l,m) */
      for (l = m + 2; l <= lmax; ++l)
        {
          k += l;
          plm = ((2*l - 1) * x * pm1 - (l + m - 1) * pm2) /
                (double) (l - m);
          result_array[k] = plm;
          result_deriv_array[k] = ((l + m) * pm1 - l * x * plm) / u_sq;
          pm2 = pm1;
          pm1 = plm;
        }
    } /* for (m = 1; m <= lmax - 1; ++m) */

  /* compute P(lmax,lmax) */

  kstart += m + 1;
  twomm1 += 2.0;
  pmm *= w->csphase * u * twomm1;
  result_array[kstart] = pmm;
  result_deriv_array[kstart] = - lmax * x * pmm / u_sq;

  return 0;
} /* gsl_Plm_deriv_array_none() */
#endif

/*
gsl_sf_legendre_array_schmidt_e()
  This routine computes Schmidt semi-normalized associated
Legendre polynomials

Inputs: lmax         - maximum order
        x            - legendre argument in [-1,1]
        csphase      - -1.0 to include CS phase (-1)^m, 1.0 to not include
        result_array - (output) where to store P_{lm}(x) values

Notes:
1) The end of the array result_array is used to store square root factors
needed in the recurrence; the front of the array stores the final
S_{lm} values
*/

int
gsl_sf_legendre_array_schmidt_e(const size_t lmax, const double x,
                                const double csphase, double result_array[])
{
  if (x > 1.0 || x < -1.0)
    {
      GSL_ERROR("x is outside [-1,1]", GSL_EDOM);
    }
  else if (csphase != 1.0 && csphase != -1.0)
    {
      GSL_ERROR("csphase has invalid value", GSL_EDOM);
    }
  else
    {
      const double eps = 1.0e-280;
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      size_t l, m;
      size_t k, idxmm;
      double plm, /* eps * S(l,m) / u^m */
             pmm; /* eps * S(m,m) / u^m */
      double rescalem;
      double pm1, /* S(l-1,m) */
             pm2; /* S(l-2,m) */
      size_t nlm = gsl_sf_legendre_nlm(lmax);
      double *sqrts = &(result_array[nlm]);

      /* precompute square root factors for recurrence */
      legendre_sqrts(lmax, sqrts);

      /* initial values S(0,0) and S(1,0) */

      pm2 = 1.0; /* S(0,0) */
      pm1 = x;   /* S(1,0) */

      result_array[0] = pm2;

      if (lmax == 0)
        return GSL_SUCCESS;

      result_array[1] = pm1;

      /* Compute S(l,0) for l=2..lmax, no scaling required */

      k = 1; /* idx(1,0) */
      for (l = 2; l <= lmax; ++l)
        {
          double linv = 1.0 / (double)l;

          k += l;  /* idx(l,m) = idx(l-1,m) + l */

          plm = (2.0 - linv) * x * pm1 - (1.0 - linv) * pm2;
          result_array[k] = plm;
          pm2 = pm1;
          pm1 = plm;
        }

      /* Compute S(m,m), S(m+1,m) and S(l,m) */

      /*
       * pi_m = Prod_{i=2}^m sqrt[ (2m - 1) / (2m) ]
       * but pi_1 = 1.0, so initialize to sqrt(2) so that
       * the first m = 1 iteration of the loop will reset it
       * to 1.0. Starting with m = 2 it will begin accumulating
       * the correct terms.
       *
       * pmm = S(m,m) * eps / u^m = pi_m
       */
      pmm = sqrt(2.0) * eps;

      rescalem = 1.0 / eps;
      idxmm = 0; /* tracks idx(m,m), initialize to idx(0,0) */

      for (m = 1; m < lmax; ++m)
        {
          /* rescalem = u^m / eps */
          rescalem *= u;

          /* compute: S(m,m) = u * sqrt((2m - 1) / (2m)) S(m-1,m-1) = u^m * pi_m */

          idxmm += m + 1; /* idx(m,m) = idx(m-1,m-1) + m + 1 */
          pmm *= csphase * sqrts[2 * m - 1] / sqrts[2 * m]; /* S(m,m) * eps / u^m */
          result_array[idxmm] = pmm * rescalem;
          pm2 = pmm;

          /* compute: S(m+1,m) = sqrt(2 * m + 1) * x * S(m,m) */

          k = idxmm + m + 1; /* idx(m+1,m) = idx(m,m) + m + 1 */
          pm1 = x * sqrts[2 * m + 1] * pm2;
          result_array[k] = pm1 * rescalem;

          /* compute S(l,m) for l=m+2...lmax */
          for (l = m + 2; l <= lmax; ++l)
            {
              k += l; /* idx(l,m) = idx(l-1,m) + l */
              plm =
                (2*l - 1) / sqrts[l + m] / sqrts[l - m] * x * pm1 -
                sqrts[l - m - 1] * sqrts[l + m - 1] /
                sqrts[l + m] / sqrts[l - m] * pm2;
              result_array[k] = plm * rescalem;
              pm2 = pm1;
              pm1 = plm;
            }
        } /* for (m = 1; m < lmax; ++m) */

      /* compute S(lmax,lmax) */

      rescalem *= u;
      idxmm += m + 1; /* idx(lmax,lmax) */
      pmm *= csphase * sqrts[2 * lmax - 1] / sqrts[2 * lmax];
      result_array[idxmm] = pmm * rescalem;

      return GSL_SUCCESS;
    }
} /* gsl_sf_legendre_array_schmidt_e() */

/*
gsl_sf_legendre_deriv_alt_array_schmidt_e()
  This routine computes Schmidt semi-normalized associated
Legendre polynomials and their first derivatives. First derivatives
of the form

d/dtheta P_{lm}(x)

are computed rather than P'_{lm}(x).

Inputs: lmax               - maximum order
        x                  - legendre argument in [-1,1]
        csphase            - -1.0 to include CS phase (-1)^m, 1.0 to not include
        result_array       - (output) where to store P_{lm}(x) values
        result_deriv_array - (output) where to store d/dtheta P_{lm}(x) values
*/

int
gsl_sf_legendre_deriv_alt_array_schmidt_e(const size_t lmax, const double x,
                                          const double csphase, double result_array[],
                                          double result_deriv_array[])
{
  if (x > 1.0 || x < -1.0)
    {
      GSL_ERROR("x is outside [-1,1]", GSL_EDOM);
    }
  else if (fabs(x) == 1.0)
    {
      GSL_ERROR("x cannot equal 1 or -1 for derivative computation", GSL_EDOM);
    }
  else if (csphase != 1.0 && csphase != -1.0)
    {
      GSL_ERROR("csphase has invalid value", GSL_EDOM);
    }
  else
    {
      const double eps = 1.0e-280;
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const double uinv = 1.0 / u;
      const double xbyu = x * uinv; /* x / u */
      size_t l, m;
      size_t k, idxmm;
      double plm, /* eps * S(l,m) / u^m */
             pmm; /* eps * S(m,m) / u^m */
      double rescalem;
      double pm1, /* S(l-1,m) */
             pm2; /* S(l-2,m) */
      size_t nlm = gsl_sf_legendre_nlm(lmax);
      double *sqrts = &(result_array[nlm]);

      /* precompute square root factors for recurrence */
      legendre_sqrts(lmax, sqrts);

      /* initial values S(0,0) and S(1,0) */
      pm2 = 1.0; /* S(0,0) */
      pm1 = x;   /* S(1,0) */

      result_array[0] = pm2;
      result_deriv_array[0] = 0.0;

      if (lmax == 0)
        return GSL_SUCCESS;

      result_array[1] = pm1;
      result_deriv_array[1] = -u;

      /* Compute S(l,0) for l=2..lmax, no scaling required */

      k = 1; /* idx(1,0) */
      for (l = 2; l <= lmax; ++l)
        {
          double linv = 1.0 / (double)l;

          k += l;  /* idx(l,m) = idx(l-1,m) + l */

          plm = (2.0 - linv) * x * pm1 - (1.0 - linv) * pm2;
          result_array[k] = plm;
          result_deriv_array[k] = uinv * l * (x * plm - pm1);
          pm2 = pm1;
          pm1 = plm;
        }

      /* Compute S(m,m), S(m+1,m) and S(l,m) */

      /*
       * pi_m = Prod_{i=2}^m sqrt[ (2m - 1) / (2m) ]
       * but pi_1 = 1.0, so initialize to sqrt(2) so that
       * the first m = 1 iteration of the loop will reset it
       * to 1.0. Starting with m = 2 it will begin accumulating
       * the correct terms.
       *
       * pmm = S(m,m) * eps / u^m = pi_m
       */
      pmm = sqrt(2.0) * eps;

      rescalem = 1.0 / eps;
      idxmm = 0; /* tracks idx(m,m), initialize to idx(0,0) */

      for (m = 1; m < lmax; ++m)
        {
          /* rescalem = u^m / eps */
          rescalem *= u;

          /*
           * compute:
           * S(m,m) = u * sqrt((2m - 1) / (2m)) S(m-1,m-1) = u^m * pi_m
           * d_t S(m,m) = m * x / u * S(m,m)
           */

          idxmm += m + 1; /* idx(m,m) = idx(m-1,m-1) + m + 1 */
          pmm *= csphase * sqrts[2 * m - 1] / sqrts[2 * m]; /* S(m,m) * eps / u^m */
          result_array[idxmm] = pmm * rescalem;
          result_deriv_array[idxmm] = m * xbyu * result_array[idxmm];
          pm2 = pmm;

          /*
           * compute:
           * S(m+1,m) = sqrt(2 * m + 1) * x * S(m,m)
           * d_t S(m+1,m) = 1/u * ((m+1)*x*S(m+1,m) - sqrt(2*m+1)*S(m,m))
           */

          k = idxmm + m + 1; /* idx(m+1,m) = idx(m,m) + m + 1 */
          pm1 = x * sqrts[2 * m + 1] * pm2;
          result_array[k] = pm1 * rescalem;
          result_deriv_array[k] =
            uinv * ((m + 1.0) * x * result_array[k] -
                    sqrts[2 * m + 1] * result_array[idxmm]);

          /* compute S(l,m) for l=m+2...lmax */
          for (l = m + 2; l <= lmax; ++l)
            {
              k += l; /* idx(l,m) = idx(l-1,m) + l */
              plm =
                (2*l - 1) / sqrts[l + m] / sqrts[l - m] * x * pm1 -
                sqrts[l - m - 1] * sqrts[l + m - 1] /
                sqrts[l + m] / sqrts[l - m] * pm2;
              result_array[k] = plm * rescalem;
              result_deriv_array[k] =
                uinv * (l * x * result_array[k] -
                        sqrts[l + m] * sqrts[l - m] * result_array[k - l]);
              pm2 = pm1;
              pm1 = plm;
            }
        } /* for (m = 1; m < lmax; ++m) */

      /* compute S(lmax,lmax) */

      rescalem *= u;
      idxmm += m + 1; /* idx(lmax,lmax) */
      pmm *= csphase * sqrts[2 * lmax - 1] / sqrts[2 * lmax];
      result_array[idxmm] = pmm * rescalem;
      result_deriv_array[idxmm] = lmax * xbyu * result_array[idxmm];

      return GSL_SUCCESS;
    }
} /* gsl_sf_legendre_deriv_alt_array_schmidt_e() */

/*
gsl_sf_legendre_deriv2_alt_array_schmidt_e()
  This routine computes Schmidt semi-normalized associated
Legendre polynomials and their first and second derivatives.

Inputs: lmax                - maximum order
        x                   - legendre argument in [-1,1]
        csphase             - -1.0 to include CS phase (-1)^m, 1.0 to not include
        result_array        - (output) where to store P_{lm}(x) values
        result_deriv_array  - (output) where to store d/dtheta P_{lm}(x) values
        result_deriv2_array - (output) where to store d^2/dtheta^2 P_{lm}(x) values
*/

int
gsl_sf_legendre_deriv2_alt_array_schmidt_e(const size_t lmax, const double x,
                                          const double csphase, double result_array[],
                                          double result_deriv_array[],
                                          double result_deriv2_array[])
{
  if (x > 1.0 || x < -1.0)
    {
      GSL_ERROR("x is outside [-1,1]", GSL_EDOM);
    }
  else if (fabs(x) == 1.0)
    {
      GSL_ERROR("x cannot equal 1 or -1 for derivative computation", GSL_EDOM);
    }
  else if (csphase != 1.0 && csphase != -1.0)
    {
      GSL_ERROR("csphase has invalid value", GSL_EDOM);
    }
  else
    {
      const double eps = 1.0e-280;
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const double uinv = 1.0 / u;
      const double uinv2 = uinv * uinv;
      const double xbyu = x * uinv; /* x / u */
      size_t l, m;
      size_t k, idxmm;
      double plm, /* eps * S(l,m) / u^m */
             pmm; /* eps * S(m,m) / u^m */
      double rescalem;
      double pm1, /* S(l-1,m) */
             pm2; /* S(l-2,m) */
      size_t nlm = gsl_sf_legendre_nlm(lmax);
      double *sqrts = &(result_array[nlm]);

      /* precompute square root factors for recurrence */
      legendre_sqrts(lmax, sqrts);

      /* initial values S(0,0) and S(1,0) */
      pm2 = 1.0; /* S(0,0) */
      pm1 = x;   /* S(1,0) */

      result_array[0] = pm2;
      result_deriv_array[0] = 0.0;
      result_deriv2_array[0] = 0.0;

      if (lmax == 0)
        return GSL_SUCCESS;

      result_array[1] = pm1;
      result_deriv_array[1] = -u;
      result_deriv2_array[1] = -x;

      /* Compute S(l,0) for l=2..lmax, no scaling required */

      k = 1; /* idx(1,0) */
      for (l = 2; l <= lmax; ++l)
        {
          double linv = 1.0 / (double)l;

          k += l;  /* idx(l,m) = idx(l-1,m) + l */

          plm = (2.0 - linv) * x * pm1 - (1.0 - linv) * pm2;
          result_array[k] = plm;
          result_deriv_array[k] = uinv * l * (x * plm - pm1);
          result_deriv2_array[k] = -(double) l * (l + 1.0) * plm -
                                   xbyu * result_deriv_array[k];
          pm2 = pm1;
          pm1 = plm;
        }

      /* Compute S(m,m), S(m+1,m) and S(l,m) */

      /*
       * pi_m = Prod_{i=2}^m sqrt[ (2m - 1) / (2m) ]
       * but pi_1 = 1.0, so initialize to sqrt(2) so that
       * the first m = 1 iteration of the loop will reset it
       * to 1.0. Starting with m = 2 it will begin accumulating
       * the correct terms.
       *
       * pmm = S(m,m) * eps / u^m = pi_m
       */
      pmm = sqrt(2.0) * eps;

      rescalem = 1.0 / eps;
      idxmm = 0; /* tracks idx(m,m), initialize to idx(0,0) */

      for (m = 1; m < lmax; ++m)
        {
          /* rescalem = u^m / eps */
          rescalem *= u;

          /*
           * compute:
           * S(m,m) = u * sqrt((2m - 1) / (2m)) S(m-1,m-1) = u^m * pi_m
           * d_t S(m,m) = m * x / u * S(m,m)
           */

          idxmm += m + 1; /* idx(m,m) = idx(m-1,m-1) + m + 1 */
          pmm *= csphase * sqrts[2 * m - 1] / sqrts[2 * m]; /* S(m,m) * eps / u^m */
          result_array[idxmm] = pmm * rescalem;
          result_deriv_array[idxmm] = m * xbyu * result_array[idxmm];
          result_deriv2_array[idxmm] =
            m * (uinv2 * m - (m + 1.0)) * result_array[idxmm] -
            xbyu * result_deriv_array[idxmm];
          pm2 = pmm;

          /*
           * compute:
           * S(m+1,m) = sqrt(2 * m + 1) * x * S(m,m)
           * d_t S(m+1,m) = 1/u * ((m+1)*x*S(m+1,m) - sqrt(2*m+1)*S(m,m))
           */

          k = idxmm + m + 1; /* idx(m+1,m) = idx(m,m) + m + 1 */
          pm1 = x * sqrts[2 * m + 1] * pm2;
          result_array[k] = pm1 * rescalem;
          result_deriv_array[k] =
            uinv * ((m + 1.0) * x * result_array[k] -
                    sqrts[2 * m + 1] * result_array[idxmm]);
          result_deriv2_array[k] =
            (m * m * uinv2 - (m + 1.0) * (m + 2.0)) * result_array[k] -
            xbyu * result_deriv_array[k];

          /* compute S(l,m) for l=m+2...lmax */
          for (l = m + 2; l <= lmax; ++l)
            {
              k += l; /* idx(l,m) = idx(l-1,m) + l */
              plm =
                (2*l - 1) / sqrts[l + m] / sqrts[l - m] * x * pm1 -
                sqrts[l - m - 1] * sqrts[l + m - 1] /
                sqrts[l + m] / sqrts[l - m] * pm2;
              result_array[k] = plm * rescalem;
              result_deriv_array[k] =
                uinv * (l * x * result_array[k] -
                        sqrts[l + m] * sqrts[l - m] * result_array[k - l]);
              result_deriv2_array[k] =
                (m * m * uinv2 - l * (l + 1.0)) * result_array[k] -
                xbyu * result_deriv_array[k];
              pm2 = pm1;
              pm1 = plm;
            }
        } /* for (m = 1; m < lmax; ++m) */

      /* compute S(lmax,lmax) */

      rescalem *= u;
      idxmm += m + 1; /* idx(lmax,lmax) */
      pmm *= csphase * sqrts[2 * lmax - 1] / sqrts[2 * lmax];
      result_array[idxmm] = pmm * rescalem;
      result_deriv_array[idxmm] = lmax * xbyu * result_array[idxmm];
      result_deriv2_array[idxmm] =
        lmax * (uinv2 * lmax - (lmax + 1.0)) * result_array[idxmm] -
        xbyu * result_deriv_array[idxmm];

      return GSL_SUCCESS;
    }
} /* gsl_sf_legendre_deriv2_alt_array_schmidt_e() */

/*
legendre_sqrts()
  Precompute square root factors needed for Legendre recurrence.
On output, array[i] = sqrt(i)
*/

static void
legendre_sqrts(const size_t lmax, double *array)
{
  size_t l;
  for (l = 0; l <= 2 * lmax + 1; ++l)
    array[l] = sqrt((double) l);
}
