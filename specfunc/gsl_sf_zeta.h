/* specfunc/gsl_sf_zeta.h
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

/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_ZETA_H__
#define __GSL_SF_ZETA_H__

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


/* Riemann Zeta Function
 * zeta(n) = Sum[ k^(-n), {k,1,Infinity} ]
 *
 * n=integer, n != 1
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zeta_int_impl(const int n, gsl_sf_result * result);
int gsl_sf_zeta_int_e(const int n, gsl_sf_result * result);


/* Riemann Zeta Function
 * zeta(x) = Sum[ k^(-s), {k,1,Infinity} ], s != 1.0
 *
 * s != 1.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zeta_impl(const double s, gsl_sf_result * result);
int gsl_sf_zeta_e(const double s, gsl_sf_result * result);


/* Hurwitz Zeta Function
 * zeta(s,q) = Sum[ (k+q)^(-s), {k,0,Infinity} ]
 *
 * s > 1.0, q > 0.0
 * exceptions: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_hzeta_impl(const double s, const double q, gsl_sf_result * result);
int gsl_sf_hzeta_e(const double s, const double q, gsl_sf_result * result);


/* Eta Function
 * eta(n) = (1-2^(1-n)) zeta(n)
 *
 * exceptions: GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_eta_int_impl(int n, gsl_sf_result * result);
int gsl_sf_eta_int_e(const int n, gsl_sf_result * result);


/* Eta Function
 * eta(s) = (1-2^(1-s)) zeta(s)
 *
 * exceptions: GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_eta_impl(const double s, gsl_sf_result * result);
int gsl_sf_eta_e(const double s, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_ZETA_H__ */
