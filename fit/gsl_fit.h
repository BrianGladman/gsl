/* fit/gsl_fit.h
 * 
 * Copyright (C) 2000 Brian Gough
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

#ifndef __GSL_FIT_H__
#define __GSL_FIT_H__

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

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

int gsl_fit_linear (const double * x, const size_t xstride,
                    const double * y, const size_t ystride,
                    size_t n,
                    double * c0, double * c1, 
                    double * cov00, double * cov01, double * cov11, 
                    double * sumsq);


int gsl_fit_wlinear (const double * x, const size_t xstride,
                     const double * w, const size_t wstride,
                     const double * y, const size_t ystride,
                     size_t n,
                     double * c0, double * c1, 
                     double * cov00, double * cov01, double * cov11, 
                     double * chisq);

int
gsl_fit_linear_est (double x, 
                    double c0, double c1, 
                    double c00, double c01, double c11,
                    double *y, double *y_err);


int gsl_fit_mul (const double * x, const size_t xstride,
                 const double * y, const size_t ystride,
                 size_t n,
                 double * c1, 
                 double * cov11, 
                 double * sumsq);

int gsl_fit_wmul (const double * x, const size_t xstride,
                  const double * w, const size_t wstride,
                  const double * y, const size_t ystride,
                  size_t n,
                  double * c1, 
                  double * cov11, 
                  double * sumsq);


int
gsl_fit_mul_est (double x, 
                 double c1, 
                 double c11,
                 double *y, double *y_err);


int
gsl_fit_wmultilinear (gsl_matrix * X,
                      const gsl_vector * w,
                      const gsl_vector * y,
                      gsl_vector * c,
                      gsl_matrix * cov,
                      double * chisq);


/* choose better names!! */

int gsl_fit_poly (const double * x, 
                  const double * w,
                  const double * y, 
                  size_t n,
                  double * c, size_t m,
                  double * chisq);

int gsl_fit_fns (const double * A, 
                 const double * w,
                 const double * y, 
                 size_t n,
                 double * c, size_t m,
                 double * chisq);

int gsl_fit_linear_nd (double * m, double * y, double * w);


__END_DECLS

#endif /* __GSL_FIT_H__ */
