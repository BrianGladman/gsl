/* Some things common to all the Monte-Carlo implementations */
/* Author: MJB */
/* RCS: $Id$ */

#ifndef __GSL_MONTE_H__
#define __GSL_MONTE_H__

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

/* Hide the function type in a typedef so that we can use it in all our
   integration functions, and make it easy to change things.
*/

typedef double (*gsl_monte_f_T)(double *);


__END_DECLS

#endif /* __GSL_MONTE_H__ */
