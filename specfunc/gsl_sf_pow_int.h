/* specfunc/gsl_sf_pow_int.h
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
#ifndef __GSL_SF_POW_INT_H__
#define __GSL_SF_POW_INT_H__

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


/* Calculate small integer powers.
 * Does not check for overflow/underflow.
 */
double gsl_sf_pow_2(const double x);
double gsl_sf_pow_3(const double x);
double gsl_sf_pow_4(const double x);
double gsl_sf_pow_5(const double x);
double gsl_sf_pow_6(const double x);
double gsl_sf_pow_7(const double x);
double gsl_sf_pow_8(const double x);
double gsl_sf_pow_9(const double x);


/* Calculate x^n.
 * Does not check for overflow/underflow.
 */
int     gsl_sf_pow_int_e(double x, int n, gsl_sf_result * result);
double  gsl_sf_pow_int(const double x, const int n);


#ifdef HAVE_INLINE
extern inline double gsl_sf_pow_2(const double x) { return x*x;   }
extern inline double gsl_sf_pow_3(const double x) { return x*x*x; }
extern inline double gsl_sf_pow_4(const double x) { double x2 = x*x;   return x2*x2;	}
extern inline double gsl_sf_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
extern inline double gsl_sf_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
extern inline double gsl_sf_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
extern inline double gsl_sf_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
extern inline double gsl_sf_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif /* HAVE_INLINE */


__END_DECLS

#endif /* __GSL_SF_POW_INT_H__ */
