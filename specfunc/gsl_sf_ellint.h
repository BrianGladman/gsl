/* specfunc/gsl_sf_ellint.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author: G. Jungman
 * RCS: $Id$
 */
#ifndef __GSL_SF_ELLINT_H__
#define __GSL_SF_ELLINT_H__

#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf_result.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


/* Legendre form of complete elliptic integrals
 *
 * K(k) = Integral[1/Sqrt[1 - k^2 Sin[t]^2], {t, 0, Pi/2}]
 * E(k) = Integral[  Sqrt[1 - k^2 Sin[t]^2], {t, 0, Pi/2}]
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_ellint_Kcomp_impl(double k, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_ellint_Kcomp_e(double k, gsl_mode_t mode, gsl_sf_result * result);

int gsl_sf_ellint_Ecomp_impl(double k, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_ellint_Ecomp_e(double k, gsl_mode_t mode, gsl_sf_result * result);


/* Legendre form of incomplete elliptic integrals
 *
 * F(phi,k)   = Integral[1/Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
 * E(phi,k)   = Integral[  Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
 * P(phi,k,n) = Integral[(1 + n Sin[t]^2)^(-1)/Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_ellint_F_impl(double phi, double k, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_ellint_F_e(double phi, double k, gsl_mode_t mode, gsl_sf_result * result);

int gsl_sf_ellint_E_impl(double phi, double k, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_ellint_E_e(double phi, double k, gsl_mode_t mode, gsl_sf_result * result);

int gsl_sf_ellint_P_impl(double phi, double k, double n, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_ellint_P_e(double phi, double k, double n, gsl_mode_t mode, gsl_sf_result * result);

int gsl_sf_ellint_D_impl(double phi, double k, double n, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_ellint_D_e(double phi, double k, double n, gsl_mode_t mode, gsl_sf_result * result);


/* Carlson's symmetric basis of functions
 *
 * RC(x,y)   = 1/2 Integral[(t+x)^(-1/2) (t+y)^(-1)], {t,0,Inf}]
 * RD(x,y,z) = 3/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-3/2), {t,0,Inf}]
 * RF(x,y,z) = 1/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2), {t,0,Inf}]
 * RJ(x,y,z,p) = 3/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2) (t+p)^(-1), {t,0,Inf}]
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_ellint_RC_impl(double x, double y, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_ellint_RC_e(double x, double y, gsl_mode_t mode, gsl_sf_result * result);

int gsl_sf_ellint_RD_impl(double x, double y, double z, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_ellint_RD_e(double x, double y, double z, gsl_mode_t mode, gsl_sf_result * result);

int gsl_sf_ellint_RF_impl(double x, double y, double z, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_ellint_RF_e(double x, double y, double z, gsl_mode_t mode, gsl_sf_result * result);

int gsl_sf_ellint_RJ_impl(double x, double y, double z, double p, gsl_mode_t mode, gsl_sf_result * result);
int gsl_sf_ellint_RJ_e(double x, double y, double z, double p, gsl_mode_t mode, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_ELLINT_H__ */
