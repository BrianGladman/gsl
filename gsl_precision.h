/* Author:  B. Gough and G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__


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


#endif /* __GSL_PRECISION_H__ */
