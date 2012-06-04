/*
** Copyright (C) 2004 Jonathan G. Underwood <j.underwood@open.ac.uk>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 3 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software 
** Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* Things to check: Not sure I'm using the correct error macros - should I be
   using the ones in the GSL manual? Based the ones here on what i see in
   coupling.c. */

#include <config.h>
#include <stdlib.h>
#include <error.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_wigner.h>

/* Static prototypes. */
static int istriangle (const int two_ja, const int two_jb, const int two_jc);
static int lndelta_e (const int two_j1, const int two_j2, const int two_j3, 
		      gsl_sf_result * lndelta);

/* This macro returns 1 if n is even, -1 if n is odd. */
#define PHASE(n) (GSL_IS_ODD(n) ? -1.0 : 1.0)

static int
istriangle (const int two_ja, const int two_jb, const int two_jc)
     /* Returns 0 if the triangle condition is met, and the sum of twice the
        angular momenta is even. Arguments are twice the value of the angular
        momenta. */
{
  if ((two_jc <= two_ja + two_jb) && (two_jc >= abs (two_ja - two_jb))
      && (GSL_IS_EVEN (two_ja + two_jb + two_jc)))
    return 1;
  else
    return 0;
}

static int
lndelta_e (const int two_j1, const int two_j2,
           const int two_j3, gsl_sf_result * lndelta)
     /* Calculates the natural log of the delta function - Zare Eq. A-2. Note:
        values returned from gsl_sf_ln_fact_e are always positive. */
{
  gsl_sf_result a1, a2, a3, a4;
  int status;

  status = gsl_sf_lnfact_e ((two_j1 + two_j2 - two_j3) / 2, &a1);
  status += gsl_sf_lnfact_e ((two_j1 - two_j2 + two_j3) / 2, &a2);
  status += gsl_sf_lnfact_e ((-two_j1 + two_j2 + two_j3) / 2, &a3);
  status += gsl_sf_lnfact_e ((two_j1 + two_j2 + two_j3) / 2 + 1, &a4);

  if (status != GSL_SUCCESS)
    OVERFLOW_ERROR (lndelta);

  lndelta->val = 0.5 * (a1.val + a2.val + a3.val - a4.val);
  lndelta->err = 0.5 * (a1.err + a2.err + a3.err + a4.err);
  lndelta->err += 2.0 * GSL_DBL_EPSILON * (a1.val + a2.val + a3.val + a4.val);

  return GSL_SUCCESS;
}

int
gsl_sf_wigner_3j_e (const int two_j1, const int two_j2, const int two_j3,
                    const int two_m1, const int two_m2, const int two_m3,
                    gsl_sf_result * result)
     /* Returns the value of the wigner 3j symbol - Zare Eq. A-1. Note that in
        Zare's book, this equation contains typos. */
{
  int v, vmin, vmax, t1, t2, t3, t4, t5, status;
  gsl_sf_result a, a1, a2, a3, a4, a5, a6, a7, a8;

  if ((two_j1 < 0) || (two_j2 < 0) || (two_j3 < 0) ||
      (abs (two_m1) > two_j1) || (abs (two_m2) > two_j2) ||
      (abs (two_m3) > two_j3) || GSL_IS_ODD (two_j1 + two_m1) ||
      GSL_IS_ODD (two_j2 + two_m2) || GSL_IS_ODD (two_j3 + two_m3))
    DOMAIN_ERROR (result);

  result->val = 0.0;
  result->err = 0.0;

  if ((!istriangle (two_j1, two_j2, two_j3)) ||
      (two_m1 + two_m2 + two_m3) != 0)
    return GSL_SUCCESS;

  t1 = (two_j1 + two_j2 - two_j3) / 2;
  t2 = (two_j1 - two_m1) / 2;
  t3 = (two_j2 + two_m2) / 2;
  t4 = (two_j3 - two_j2 + two_m1) / 2;
  t5 = (two_j3 - two_j1 - two_m2) / 2;

  vmin = (t4 < t5) ? -t4 : -t5;
  vmin = (vmin < 0) ? 0 : vmin;

  vmax = t1;
  vmax = (t2 < vmax) ? t2 : vmax;
  vmax = (t3 < vmax) ? t3 : vmax;

  if (vmin > vmax)
    return GSL_SUCCESS;

  status = gsl_sf_lnfact_e ((two_j1 + two_m1) / 2, &a1);
  status += gsl_sf_lnfact_e (t2, &a2);
  status += gsl_sf_lnfact_e (t3, &a3);
  status += gsl_sf_lnfact_e ((two_j2 - two_m2) / 2, &a4);
  status += gsl_sf_lnfact_e ((two_j3 + two_m3) / 2, &a5);
  status += gsl_sf_lnfact_e ((two_j3 - two_m3) / 2, &a6);
  status += lndelta_e (two_j1, two_j2, two_j3, &a7);

  if (status != GSL_SUCCESS)
    OVERFLOW_ERROR (result);

  a.val =
    0.5 * (a1.val + a2.val + a3.val + a4.val + a5.val + a6.val) + a7.val;
  a.err =
    0.5 * (a1.err + a2.err + a3.err + a4.err + a5.err + a6.err) + a7.err;
  a.err += 2.0 * GSL_DBL_EPSILON * a.val;

  for (v = vmin; v <= vmax; v++)
    {
      status += gsl_sf_lnfact_e (v, &a1);
      status += gsl_sf_lnfact_e (t1 - v, &a2);
      status += gsl_sf_lnfact_e (t2 - v, &a3);
      status += gsl_sf_lnfact_e (t3 - v, &a4);
      status += gsl_sf_lnfact_e (t4 + v, &a5);
      status += gsl_sf_lnfact_e (t5 + v, &a6);
      if (status != GSL_SUCCESS)
        OVERFLOW_ERROR (result);

      a7.val = a.val - (a1.val + a2.val + a3.val + a4.val + a5.val + a6.val);
      a7.err = (a.err + a1.err + a2.err + a3.err + a4.err + a5.err + a6.err);
      a7.err += 2.0 * GSL_DBL_EPSILON *
        (a.val + a1.val + a2.val + a3.val + a4.val + a5.val + a6.val);

      status += gsl_sf_exp_err_e (a7.val, a7.err, &a8);
      if (status != GSL_SUCCESS)
        OVERFLOW_ERROR (result);

      result->val += PHASE (v) * a8.val;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs (result->val);
    }

  result->val *= PHASE ((two_j1 - two_j2 - two_m3) / 2);
  result->err += 2.0 * GSL_DBL_EPSILON * fabs (result->val);

  return GSL_SUCCESS;
}

int
gsl_sf_wigner_6j_e (const int two_j1, const int two_j2, const int two_j3,
                    const int two_j4, const int two_j5, const int two_j6,
                    gsl_sf_result * result)
     /* Returns the value of the 6j symbol - Zare Eq. A-4 */
{
  int v, vmin, vmax, t1, t2, t3, t4, t5, t6, t7, status;
  gsl_sf_result a, a1, a2, a3, a4;

  if ((two_j1 < 0) || (two_j2 < 0) || (two_j3 < 0) ||
      (two_j4 < 0) || (two_j5 < 0) || (two_j6 < 0))
    DOMAIN_ERROR (result);

  result->val = 0.0;
  result->err = 0.0;

  if (!istriangle (two_j1, two_j2, two_j3) ||
      !istriangle (two_j1, two_j5, two_j6) ||
      !istriangle (two_j4, two_j2, two_j6) ||
      !istriangle (two_j4, two_j5, two_j3))
    return GSL_SUCCESS;

  t1 = (two_j1 + two_j2 + two_j3) / 2;
  t2 = (two_j1 + two_j5 + two_j6) / 2;
  t3 = (two_j4 + two_j2 + two_j6) / 2;
  t4 = (two_j4 + two_j5 + two_j3) / 2;
  t5 = (two_j1 + two_j2 + two_j4 + two_j5) / 2;
  t6 = (two_j2 + two_j3 + two_j5 + two_j6) / 2;
  t7 = (two_j3 + two_j1 + two_j6 + two_j4) / 2;

  vmin = t1;
  vmin = (t2 > vmin) ? t2 : vmin;
  vmin = (t3 > vmin) ? t3 : vmin;
  vmin = (t4 > vmin) ? t4 : vmin;

  vmax = t5;
  vmax = (t6 < vmax) ? t6 : vmax;
  vmax = (t7 < vmax) ? t7 : vmax;

  if (vmin > vmax)
    return GSL_SUCCESS;

  status = lndelta_e (two_j1, two_j2, two_j3, &a1);
  status += lndelta_e (two_j1, two_j5, two_j6, &a2);
  status += lndelta_e (two_j4, two_j2, two_j6, &a3);
  status += lndelta_e (two_j4, two_j5, two_j3, &a4);

  if (status != GSL_SUCCESS)
    OVERFLOW_ERROR (result);

  a.val = a1.val + a2.val + a3.val + a4.val;
  a.err = a1.err + a2.err + a3.err + a4.err;
  a.err += 2.0 * GSL_DBL_EPSILON * a.val;

  for (v = vmin; v <= vmax; v++)
    {
      gsl_sf_result b1, b2, b3, b4, b5, b6, b7, b8, b9, b10;
      status += gsl_sf_lnfact_e (v + 1, &b1);
      status += gsl_sf_lnfact_e (v - t1, &b2);
      status += gsl_sf_lnfact_e (v - t2, &b3);
      status += gsl_sf_lnfact_e (v - t3, &b4);
      status += gsl_sf_lnfact_e (v - t4, &b5);
      status += gsl_sf_lnfact_e (t5 - v, &b6);
      status += gsl_sf_lnfact_e (t6 - v, &b7);
      status += gsl_sf_lnfact_e (t7 - v, &b8);
      if (status != GSL_SUCCESS)
        OVERFLOW_ERROR (result);

      b9.val = a.val + b1.val - b2.val - b3.val - b4.val - b5.val - b6.val -
        b7.val - b8.val;
      b9.err = a.err + b1.err + b2.err + b3.err + b4.err + b5.err + b6.err +
        b7.err + b8.err;
      b9.err += 2.0 * GSL_DBL_EPSILON * (a.val + b1.val + b2.val + b3.val +
                                         b4.val + b5.val + b6.val + b7.val +
                                         b8.val);

      status += gsl_sf_exp_err_e (b9.val, b9.err, &b10);
      if (status != GSL_SUCCESS)
        OVERFLOW_ERROR (result);

      result->val += b10.val * PHASE (v);
      result->err += b10.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs (result->val);
    }

  return GSL_SUCCESS;
}

int
gsl_sf_wigner_9j_e (const int two_j1, const int two_j2, const int two_j3,
                    const int two_j4, const int two_j5, const int two_j6,
                    const int two_j7, const int two_j8, const int two_j9,
                    gsl_sf_result * result)
    /* Returns the value of 9j symbol - Zare Eq. 4.24 */
{
  int two_k, two_kmin, two_kmax, status = GSL_SUCCESS;

  if ((two_j1 < 0) || (two_j2 < 0) || (two_j3 < 0) ||
      (two_j4 < 0) || (two_j5 < 0) || (two_j6 < 0) ||
      (two_j7 < 0) || (two_j8 < 0) || (two_j9 < 0))
    DOMAIN_ERROR (result);

  result->val = 0.0;
  result->err = 0.0;

  if (!istriangle (two_j1, two_j2, two_j3) ||
      !istriangle (two_j4, two_j5, two_j6) ||
      !istriangle (two_j7, two_j8, two_j9) ||
      !istriangle (two_j1, two_j4, two_j7) ||
      !istriangle (two_j2, two_j5, two_j8) ||
      !istriangle (two_j3, two_j6, two_j9))
    return GSL_SUCCESS;

  two_kmin = abs (two_j1 - two_j9);
  two_kmin =
    (abs (two_j8 - two_j4) < two_kmin) ? abs (two_j8 - two_j4) : two_kmin;
  two_kmin =
    (abs (two_j6 - two_j2) < two_kmin) ? abs (two_j6 - two_j2) : two_kmin;

  two_kmax = (two_j1 + two_j9);
  two_kmax = (two_j8 + two_j4 > two_kmax) ? two_j8 + two_j4 : two_kmax;
  two_kmax = (two_j6 + two_j2 > two_kmax) ? two_j6 + two_j2 : two_kmax;

  if (two_kmax < two_kmin)
    {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  for (two_k = two_kmin; two_k <= two_kmax; two_k += 2)
    {
      gsl_sf_result a1, a2, a3;
      status += gsl_sf_wigner_6j_e (two_j1, two_j4, two_j7,
                                    two_j8, two_j9, two_k, &a1);
      status +=
        gsl_sf_wigner_6j_e (two_j2, two_j5, two_j8,
                            two_j4, two_k, two_j6, &a2);
      status +=
        gsl_sf_wigner_6j_e (two_j3, two_j6, two_j9,
                            two_k, two_j1, two_j2, &a3);
      if (status != GSL_SUCCESS)
        OVERFLOW_ERROR (result);

      result->val += (two_k + 1.0) * PHASE (two_k) * a1.val * a2.val * a3.val;

      result->err += result->val * (fabs (a1.err / a1.val) +
                                    fabs (a2.err / a2.val) +
                                    fabs (a3.err / a3.val));
      result->err += 2.0 * GSL_DBL_EPSILON * fabs (result->val);
    }
  return GSL_SUCCESS;
}

int
gsl_sf_wigner_drot_e (const int two_j, const int two_m1, const int two_m2,
                      const double theta, gsl_sf_result * result)
     /* Returns the value of the reduced rotation matrix - Zare Eq. 3.57. (which
        contains typos).  Angle theta is in radians. */
{
  int v, vmin, vmax, t1, t2, t3, t4, tv, status = GSL_SUCCESS;
  gsl_sf_result a, a1, a2, a3, a4, st, ct;

  if ((two_j) < 0 || abs (two_m1) > two_j || abs (two_m2) > two_j ||
      (GSL_IS_ODD (two_j + two_m1)) || (GSL_IS_ODD (two_j + two_m2)))
    DOMAIN_ERROR (result);

  result->val = 0.0;
  result->err = 0.0;

  t1 = (two_j + two_m1) / 2;
  t2 = (two_j - two_m2) / 2;
  t3 = (two_m2 - two_m1) / 2;
  t4 = (2 * two_j + two_m1 - two_m2) / 2;

  vmin = (-t3 < 0) ? 0 : -t3;
  vmax = (t2 < t1) ? t2 : t1;

  if (vmin > vmax)
    return GSL_SUCCESS;

  status += gsl_sf_cos_e (theta * 0.5, &ct);
  status += gsl_sf_sin_e (theta * 0.5, &st);
  if (status != GSL_SUCCESS)
    OVERFLOW_ERROR (result);

  status += gsl_sf_lnfact_e (t1, &a1);
  status += gsl_sf_lnfact_e (t2, &a2);
  status += gsl_sf_lnfact_e ((two_j - two_m1) / 2, &a3);
  status += gsl_sf_lnfact_e ((two_j + two_m2) / 2, &a4);
  if (status != GSL_SUCCESS)
    OVERFLOW_ERROR (result);

  a.val = 0.5 * (a1.val + a2.val + a3.val + a4.val);
  a.err = 0.5 * (a1.err + a2.err + a3.err + a4.err);
  a.err += 2.0 * GSL_DBL_EPSILON * a.val;

  for (v = vmin; v <= vmax; v++)
    {
      gsl_sf_result b1, b2, b3, b4, b5, b6, b7, b8;
      int i1, i2;

      status += gsl_sf_lnfact_e (t1 - v, &b1);
      status += gsl_sf_lnfact_e (t2 - v, &b2);
      status += gsl_sf_lnfact_e (t3 + v, &b3);
      status += gsl_sf_lnfact_e (v, &b4);
      if (status != GSL_SUCCESS)
        OVERFLOW_ERROR (result);

      b5.val = a.val - b1.val - b2.val - b3.val - b4.val;
      b5.err = a.err + b1.err + b2.err + b3.err + b4.err;
      b5.err +=
        2.0 * GSL_DBL_EPSILON * (a.val + b1.val + b2.val + b3.val + b4.val);

      status += gsl_sf_exp_err_e (b5.val, b5.err, &b6);
      if (status != GSL_SUCCESS)
        OVERFLOW_ERROR (result);

      tv = 2 * v;
      i1 = t4 - tv;
      i2 = t3 + tv;

      status += gsl_sf_pow_int_e (ct.val, i1, &b7);
      status += gsl_sf_pow_int_e (st.val, i2, &b8);
      if (status != GSL_SUCCESS)
        OVERFLOW_ERROR (result);

      /* Bolt on error in pow_int for the error in the input values */
      b7.err += i1 * fabs (ct.err);
      b8.err += i2 * fabs (st.err);

      result->val += PHASE (v) * b6.val * b7.val * b8.val;

      result->err += fabs (result->val) *
        (fabs (b6.err / b6.val) + fabs (b7.err / b7.val) +
         fabs (b8.err / b8.val));
      result->err += 2.0 * GSL_DBL_EPSILON * fabs (result->val);
    }

  return GSL_SUCCESS;
}

#include "eval.h"

double
gsl_sf_wigner_3j (const int two_j1, const int two_j2, const int two_j3,
                  const int two_m1, const int two_m2, const int two_m3)
{
  EVAL_RESULT (gsl_sf_wigner_3j_e (two_j1, two_j2, two_j3,
                                   two_m1, two_m2, two_m3, &result));
}

double
gsl_sf_wigner_6j (const int two_j1, const int two_j2, const int two_j3,
                  const int two_j4, const int two_j5, const int two_j6)
{
  EVAL_RESULT (gsl_sf_wigner_6j_e (two_j1, two_j2, two_j3,
                                   two_j4, two_j5, two_j6, &result));
}

double
gsl_sf_wigner_9j (const int two_j1, const int two_j2, const int two_j3,
                  const int two_j4, const int two_j5, const int two_j6,
                  const int two_j7, const int two_j8, const int two_j9)
{
  EVAL_RESULT (gsl_sf_wigner_9j_e (two_j1, two_j2, two_j3,
                                   two_j4, two_j5, two_j6,
                                   two_j7, two_j8, two_j9, &result));
}

double
gsl_sf_wigner_drot (const int two_j, const int two_m1, const int two_m2,
                    const double theta)
{
  EVAL_RESULT (gsl_sf_wigner_drot_e (two_j, two_m1, two_m2, theta, &result));
}

#undef PHASE
