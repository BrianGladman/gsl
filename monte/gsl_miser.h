/* Author: MJB */
/* RCS: $Id$ */

extern size_t min_calls;
extern size_t min_calls_per_bisection;
  
extern double dither;

int gsl_monte_miser(double (*func)(double []), double xl[], double xh[], 
		    size_t num_dim, size_t calls, double *ave, double *var);
