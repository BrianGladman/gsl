/* Author:  B. Gough and G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__

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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
extern const double gsl_prec_eps[];
extern const double gsl_sqrt_prec_eps[];
extern const double gsl_root3_prec_eps[];
extern const double gsl_root4_prec_eps[];
extern const double gsl_root5_prec_eps[];
extern const double gsl_root6_prec_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */
