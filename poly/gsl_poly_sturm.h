/* poly/gsl_poly.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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

#ifndef __GSL_POLY_H__
#define __GSL_POLY_H__

#include <stdlib.h>
#include <gsl/gsl_complex.h>

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


/* Struct for a Sturm sequence */
struct gsl_poly_sturm_struct {

   gsl_poly ** sturmseq;
   size_t index;
   size_t size;
};
typedef struct gsl_poly_sturm_struct gsl_poly_sturm;

/* Allocate Sturm sequence */
gsl_poly_sturm *
gsl_poly_sturm_alloc(const size_t size);

/* Allocate Sturm sequence and initialize all data to zero*/
gsl_poly_sturm *
gsl_poly_sturm_calloc(const size_t size);

/* Free Sturm sequence */
void gsl_poly_sturm_free(gsl_poly_sturm * ss);


/* Build sturm sequence */
int gsl_poly_sturm_build(gsl_poly_sturm *ss, const gsl_poly * p, gsl_poly * w);

/* Calculate the number of sign changes in the Sturm sequence at point a */

int gsl_poly_sturm_changes(const gsl_poly_sturm * ss, double a);

/* Calculate the number of roots between a and b based on Sturm sequence ss */
int gsl_poly_sturm_numroots(const gsl_poly_sturm * ss,double a, double b);


/* Dump Sturm sequence */

void gsl_poly_sturm_dump(const gsl_poly_sturm * ss);




