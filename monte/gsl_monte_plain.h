/* monte/gsl_monte_plain.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
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

/* Plain Monte-Carlo. */

/* Author: MJB */
/* RCS: $Id$ */

#ifndef __GSL_MONTE_PLAIN_H__
#define __GSL_MONTE_PLAIN_H__

#include <stdio.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>

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

typedef struct {
  /* control variables */
  double acc;

  int init_done;
  int check_done;
  int verbose;

  size_t num_dim;

  FILE* ostream;
  gsl_rng* ranf;

} gsl_monte_plain_state;

int gsl_monte_plain_integrate(gsl_monte_plain_state *state, 
			      const gsl_monte_f_T fun, 
			      const double* xl, const double* xu, 
			      const size_t num_dim, 
			      const size_t calls, double* res, double* err);

gsl_monte_plain_state* gsl_monte_plain_alloc(size_t num_dim);

int gsl_monte_plain_validate(gsl_monte_plain_state* state,
                             const double xl[], const double xu[], 
                             unsigned long num_dim, unsigned long calls);

int gsl_monte_plain_init(gsl_monte_plain_state* state);

void gsl_monte_plain_free (gsl_monte_plain_state* s);

__END_DECLS

#endif /* __GSL_MONTE_PLAIN_H__ */
