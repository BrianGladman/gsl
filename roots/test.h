/* roots/test.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Reid Priedhorsky, Brian Gough
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


/* $Id$ */

typedef double simple_function (double x);
typedef struct { simple_function * f; simple_function * df; } function_pair ;

gsl_function create_function (simple_function * f) ;
double eval_function (double x, void * params) ;

gsl_function_fdf create_fdf (simple_function * f, simple_function * df);
double eval_fdf_f (double x, void * params);
double eval_fdf_df (double x, void * params);
void eval_fdf (double x, void * params, double * y1, double * y2);

void
  test_macros (void);

void
  test_roots (void);

void
  test_poly (void);

void
test_f (const gsl_root_fsolver_type * T, 
	const char * description, gsl_function *f,
	double lower_bound, double upper_bound, double correct_root);

void
test_f_e (const gsl_root_fsolver_type * T, const char * description, 
	  gsl_function *f,
	  double lower_bound, double upper_bound, double correct_root);

void
test_fdf (const gsl_root_fdfsolver_type * T, const char * description, 
	  gsl_function_fdf *fdf, double root, double correct_root);

void
test_fdf_e (const gsl_root_fdfsolver_type * T, const char * description, 
	    gsl_function_fdf *fdf, double root, double correct_root);


void
  usage (void);

void
  error_handler (const char *reason, const char *file, int line);

double
  func1 (double x);

double
  func1_df (double x);

void
  func1_fdf (double x, double *y, double *yprime);

double
  func2 (double x);

double
  func2_df (double x);

void
  func2_fdf (double x, double *y, double *yprime);

double
  func3 (double x);

double
  func3_df (double x);

void
  func3_fdf (double x, double *y, double *yprime);

double
  func4 (double x);

double
  func4_df (double x);

void
  func4_fdf (double x, double *y, double *yprime);

double
  func5 (double x);

double
  func5_df (double x);

void
  func5_fdf (double x, double *y, double *yprime);

double
  func6 (double x);

double
  func6_df (double x);

void
  func6_fdf (double x, double *y, double *yprime);

double
  sin_f (double x);

double
  sin_df (double x);

void
  sin_fdf (double x, double *y, double *yprime);

double
  cos_f (double x);

double
  cos_df (double x);

void
  cos_fdf (double x, double *y, double *yprime);
