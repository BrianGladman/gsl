/* Plain Monte-Carlo. */

/* Author: MJB */
/* RCS: $Id$ */

#ifndef GSL_MONTE_PLAIN_H
#define GSL_MONTE_PLAIN_H

#include <gsl_monte.h>
#include <gsl_rng.h>

int gsl_monte_plain(const gsl_rng *r, const gsl_monte_f_T fun, 
		    const double* xl, const double* xh, const size_t num_dim, 
		    const size_t calls, double* avg, double* var);


#endif /* GSL_MONTE_PLAIN_H */
