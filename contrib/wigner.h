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

#ifndef __GSL_SF_WIGNER_H__
#define  __GSL_SF_WIGNER_H__ 1
#include <gsl/gsl_sf_result.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS          /* empty */
# define __END_DECLS            /* empty */
#endif

__BEGIN_DECLS
  double gsl_sf_wigner_3j (const int two_j1, const int two_j2,
                           const int two_j3, const int two_m1,
                           const int two_m2, const int two_m3);
int gsl_sf_wigner_3j_e (const int two_j1, const int two_j2, const int two_j3,
                        const int two_m1, const int two_m2, const int two_m3,
                        gsl_sf_result * result);

double gsl_sf_wigner_6j (const int two_j1, const int two_j2, const int two_j3,
                         const int two_j4, const int two_j5,
                         const int two_j6);
int gsl_sf_wigner_6j_e (const int two_j1, const int two_j2, const int two_j3,
                        const int two_j4, const int two_j5, const int two_j6,
                        gsl_sf_result * result);

double gsl_sf_wigner_9j (const int two_j1, const int two_j2, const int two_j3,
                         const int two_j4, const int two_j5, const int two_j6,
                         const int two_j7, const int two_j8,
                         const int two_j9);
int gsl_sf_wigner_9j_e (const int two_j1, const int two_j2, const int two_j3,
                        const int two_j4, const int two_j5, const int two_j6,
                        const int two_j7, const int two_j8, const int two_j9,
                        gsl_sf_result * result);

double gsl_sf_wigner_drot (const int two_j, const int two_m1,
                           const int two_m2, const double theta);
int gsl_sf_wigner_drot_e (const int two_j, const int two_m1, const int two_m2,
                          const double theta, gsl_sf_result * result);

__END_DECLS
#endif /* __GSL_SF_WIGNER_H__ */
