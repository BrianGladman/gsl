/* specfunc/gsl_sf_mathieu.h
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

#ifndef _GSL_SF_MATHIEU_H_
#define _GSL_SF_MATHIEU_H_

#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_eigen.h>

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


typedef struct 
{
    size_t size;
    size_t even_order;
    size_t odd_order;
    int extra_values;
    double *char_value;
    double *dd;
    double *ee;
    double *tt;
    double *e2;
    double *zz;
    gsl_vector *eval;
    gsl_matrix *evec;
    gsl_eigen_symmv_workspace *wmat;
} gsl_sf_mathieu_workspace;


/* Allocate computational storage space for eigenvalue solution. */
gsl_sf_mathieu_workspace *gsl_sf_mathieu_alloc(const size_t nn);
void gsl_sf_mathieu_free(gsl_sf_mathieu_workspace *workspace);

/* Compute an angular Mathieu function. */
int gsl_sf_mathieu_c(int order, double qq, double zz, gsl_sf_result *result);
int gsl_sf_mathieu_s(int order, double qq, double zz, gsl_sf_result *result);
int gsl_sf_mathieu_c_array(int nmin, int nmax, double qq, double zz,
                           gsl_sf_mathieu_workspace *work,
                           double result_array[]);
int gsl_sf_mathieu_s_array(int nmin, int nmax, double qq, double zz,
                           gsl_sf_mathieu_workspace *work,
                           double result_array[]);

/* Compute a radial Mathieu function. */
int gsl_sf_mathieu_mc_1(int order, double qq, double zz,
                        gsl_sf_result *result);
int gsl_sf_mathieu_ms_1(int order, double qq, double zz,
                        gsl_sf_result *result);


__END_DECLS

#endif /* !_GSL_SF_MATHIEU_H_ */
