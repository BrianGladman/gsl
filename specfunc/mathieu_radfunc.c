/* specfunc/mathieu_radfunc.c
 * 
 * Copyright (C) 2002 Lowell Johnson
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

/* Author:  L. Johnson */

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_mathieu.h>

#include "mathieu.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b)) 


int gsl_sf_mathieu_mc_1(int order, double qq, double zz,
                        gsl_sf_result *result)
{
  int even_odd, kk, mm, status;
  double maxerr = 1e-14, amax, pi = acos(-1.0), fn;
  double coeff[NUM_MATHIEU_COEFF], aa, fc, fj, fjp;
  double j1c, j2c, j1pc, j2pc;
  double u1, u2;


  /* Check for out of bounds parameters. */
  if (qq == 0.0)
  {
      return GSL_FAILURE;
  }

  mm = 0;
  amax = 0.0;
  fn = 0.0;
/*   rm3 = (0.0,0.0); */
/*   rm1p = 0.0; */
/*   rm3p = (0.0,0.0); */
  u1 = sqrt(qq)*exp(-1.0*zz);
  u2 = sqrt(qq)*exp(zz);
  
  even_odd = 0;
  if (order % 2 != 0)
      even_odd = 1;

  /* Compute the characteristic value. */
  status = gsl_sf_mathieu_c_charv(order, qq, &aa);
  if (status != GSL_SUCCESS)
  {
      return status;
  }
  
  /* Compute the series coefficients. */
  status = gsl_sf_mathieu_c_coeff(order, qq, aa, coeff);
  if (status != GSL_SUCCESS)
  {
      return status;
  }

  if (even_odd == 0)
  {
      for (kk=0; kk<NUM_MATHIEU_COEFF; kk++)
      {
          amax = MAX(amax, fabs(coeff[kk]));
          if (fabs(coeff[kk])/amax < maxerr)
              break;

          j1c = gsl_sf_bessel_Jn(kk, u1);
          j2c = gsl_sf_bessel_Jn(kk, u2);
              
/*          besjyf(u1, k1, j1c, y1c, h11c, h21c, j1p, y1p, h11p, h21p); */
/*          besjyf(u2, k1, j2c, y2c, h12c, h22c, j2p, y2p, h12p, h22p); */

          fc = pow(-1.0, 0.5*order+kk)*coeff[kk];
          fj = fc*j1c;
          fn += fj*j2c;
/*          rm3 += fj*h12c; */
/*          fj1 = fc*j1p*u1; */
/*          fj2 = fj*u2; */
/*          rm1p += -1.0*fj1*j2c + fj2*j2p; */
/*          rm3p += -1.0*fj1*h12c + fj2*h12p; */
      }

      fn *= sqrt(pi/2.0)/coeff[0];
/*           rm3 *= sqrt(pi/2.0)/coeff[0]; */
/*           rm1p *= sqrt(pi/2.0)/coeff[0]; */
/*           rm3p *= sqrt(pi/2.0)/coeff[0]; */
  }
  else
  {
      for (kk=0; kk<NUM_MATHIEU_COEFF; kk++)
      {
          amax = MAX(amax, fabs(coeff[kk]));
          if (fabs(coeff[kk])/amax < maxerr)
              break;

          j1c = gsl_sf_bessel_Jn(kk, u1);
          j2c = gsl_sf_bessel_Jn(kk, u2);
          j1pc = gsl_sf_bessel_Jn(kk+1, u1);
          j2pc = gsl_sf_bessel_Jn(kk+1, u2);
/*          besjyf(u1, k1, j1c, y1c, h11c, h21c, j1p, y1p, h11p, h21p); */
/*          besjyf(u1, kp1, j1pc, y1pc, h11pc, h21pc, h1pp, y1pp, */
/*                 h11pp, h21pp); */
/*          besjyf(u2, k1, j2c, y2c, h12c, h22c, j2p, y2p, h12p, h22p); */
/*          besjyf(u2, kp1, j2pc, y2pc, h12pc, h22pc, j2pp, y2pp, */
/*                 h12pp, h22pp); */
              
          fc = pow(-1.0, 0.5*(order-1)+kk)*coeff[kk];
          fj = fc*j1c;
          fjp = fc*j1pc;
          fn += fj*j2pc + fjp*j2c;
/*          rm3 += fj*h12pc + fjp*h12c; */
/*          fj1 = fc*j1p*u1; */
/*          fjp1 = fc*j1pp*u1; */
/*          fj2 = fj*u2; */
/*          fjp2 = fc*j1pc*u2; */
/*          rm1p += -1.0*fj1*j2pc + fj2*j2pp - fjp1*j2c + fjp2*j2p; */
/*          rm3p += -1.0*fj1*h12pc + fj2*h12pp - fjp1*h12c + fjp2*h12p; */
      }

      fn *= sqrt(pi/2.0)/coeff[0];
/*      rm3 *= sqrt(pi/2.0)/coeff[0]; */
/*      rm1p *= sqrt(pi/2.0)/coeff[0]; */
/*      rm3p *= sqrt(pi/2.0)/coeff[0]; */
  }

  result->val = fn;
  result->err = GSL_DBL_EPSILON*fabs(fn);
  
  return GSL_SUCCESS;
}


int gsl_sf_mathieu_ms_1(int order, double qq, double zz,
                        gsl_sf_result *result)
{
  int even_odd, kk, mm, status;
  double maxerr = 1e-14, amax, pi = acos(-1.0), fn;
  double coeff[NUM_MATHIEU_COEFF], aa, fc, fj, fjp, fjm;
  double j1c, j2c, j1mc, j2mc, j1pc, j2pc;
  double u1, u2;


  /* Check for out of bounds parameters. */
  if (qq == 0.0)
  {
      return GSL_FAILURE;
  }

  mm = 0;
  amax = 0.0;
  fn = 0.0;
/*   rm3 = (0.0,0.0); */
/*   rm1p = 0.0; */
/*   rm3p = (0.0,0.0); */
  u1 = sqrt(qq)*exp(-1.0*zz);
  u2 = sqrt(qq)*exp(zz);
  
  even_odd = 0;
  if (order % 2 != 0)
      even_odd = 1;
  
  /* Compute the characteristic value. */
  status = gsl_sf_mathieu_s_charv(order, qq, &aa);
  if (status != GSL_SUCCESS)
  {
      return status;
  }
  
  /* Compute the series coefficients. */
  status = gsl_sf_mathieu_s_coeff(order, qq, aa, coeff);
  if (status != GSL_SUCCESS)
  {
      return status;
  }

  if (even_odd == 0)
  {
      for (kk=0; kk<NUM_MATHIEU_COEFF; kk++)
      {
          amax = MAX(amax, fabs(coeff[kk]));
          if (fabs(coeff[kk])/amax < maxerr)
              break;

          j1mc = gsl_sf_bessel_Jn(kk, u1);
          j2mc = gsl_sf_bessel_Jn(kk, u2);
          j1pc = gsl_sf_bessel_Jn(kk+2, u1);
          j2pc = gsl_sf_bessel_Jn(kk+2, u2);
/*           besjyf(u1, km1, j1mc, y1mc, h11mc, h21mc, j1mp, y1mp, */
/*                  h11mp, h21mp); */
/*           besjyf(u2, km1, j2mc, y2mc, h12mc, h22mc, j2mp, y2mp, */
/*                  h12mp, h22mp); */
/*           besjyf(u1, kp1, j1pc, y1pc, h11pc, h21pc, j1pp, y1pp, */
/*                  h11pp, h21pp); */
/*           besjyf(u2, kp1, j2pc, y2pc, h12pc, h22pc, j2pp, y2pp, */
/*                  h12pp, h22pp); */
              
          fc = pow(-1.0, 0.5*order+kk+1)*coeff[kk];
          fjm = fc*j1mc;
          fjp = fc*j1pc;
          fn += fjm*j2pc - fjp*j2mc;
/*           rm3 += fjm*h12pc - fjp*h12mc; */
/*           fjm1 = fc*j1mp*u1; */
/*           fjp1 = fc*j1pp*u1; */
/*           fjm2 = fjm*u2; */
/*           fjp2 = fc*j1pc*u2; */
/*           rm1p += -1.0*fjm1*j2pc + fjm2*j2pp + fjp1*j2mc - fjp2*j2mp; */
/*           rm3p += -1.0*fjm1*h12pc + fjm2*h12pp + fjp1*h12mc - fjp2*h12mp; */
      }

      fn *= sqrt(pi/2.0)/coeff[0];
/*       rm3 *= sqrt(pi/2.0)/coeff[0]; */
/*       rm1p *= sqrt(pi/2.0)/coeff[0]; */
/*       rm3p *= sqrt(pi/2.0)/coeff[0]; */
  }
  else
  {
      for (kk=0; kk<NUM_MATHIEU_COEFF; kk++)
      {
          amax = MAX(amax, fabs(coeff[kk]));
          if (fabs(coeff[kk])/amax < maxerr)
              break;

          j1c = gsl_sf_bessel_Jn(kk, u1);
          j2c = gsl_sf_bessel_Jn(kk, u2);
          j1pc = gsl_sf_bessel_Jn(kk+1, u1);
          j2pc = gsl_sf_bessel_Jn(kk+1, u2);
/*           besjyf(u1, k1, j1c, y1c, h11c, h21c, j1p, y1p, h11p, h21p); */
/*           besjyf(u2, k1, j2c, y2c, h12c, h22c, j2p, y2p, h12p, h22p); */
/*           besjyf(u1, kp1, j1pc, y1pc, h11pc, h21pc, j1pp, y1pp, */
/*                  h11pp, h21pp); */
/*           besjyf(u2, kp1, j2pc, y2pc, h12pc, h22pc, j2pp, y2pp, */
/*                  h12pp, h22pp); */

          fc = pow(-1.0, 0.5*(order-1)+kk)*coeff[kk];
          fj = fc*j1c;
          fjp = fc*j1pc;
          fn += fj*j2pc - fjp*j2c;
/*           rm3 += fj*h12pc - fjp*h12c; */
/*           fj1 = fc*j1p*u1; */
/*           fjp1 = fc*j1pp*u1; */
/*           fj2 = fj*u2; */
/*           fjp2 = fc*j1pc*u2; */
/*           rm1p += -1.0*fj1*j2pc + fj2*j2pp + fjp1*j2c - fjp2*j2p; */
/*           rm3p += -1.0*fj1*h12pc + fj2*h12pp + fjp1*h12c - fjp2*h12p; */
      }

      fn *= sqrt(pi/2.0)/coeff[0];
/*       rm3 *= sqrt(pi/2.0)/coeff[0]; */
/*       rm1p *= sqrt(pi/2.0)/coeff[0]; */
/*       rm3p *= sqrt(pi/2.0)/coeff[0]; */
  }

  result->val = fn;
  result->err = GSL_DBL_EPSILON*fabs(fn);
  
  return GSL_SUCCESS;
}
