/* Some things common to all the Monte-Carlo implementations */
/* Author: MJB */
/* RCS: $Id$ */

#ifndef __GSL_MONTE_H__
#define __GSL_MONTE_H__

/* Hide the function type in a typedef so that we can use it in all our
   integration functions, and make it easy to change things.
*/

typedef double (*gsl_monte_f_T)(double *);


#endif /* __GSL_MONTE_H__ */
