/* Author:  B. Gough and G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_MODE_H
#define GSL_MODE_H


/* Some functions can take a mode argument. This
 * is a rough method to do things like control
 * the precision of the algorithm. This mainly
 * occurs in special functions, but we figured
 * it was ok to have a general facility.
 *
 * The mode type is 32-bit field. Most of
 * the fields are currently unused. Users
 * '|' various predefined constants to get
 * a desired mode.
 */
typedef unsigned int gsl_mode_t;


/* Here are the predefined constants.
 * Note that the precision constants
 * are special because they are used
 * to index arrays, so do not change
 * them. The precision information is
 * in the low order 3 bits of gsl_mode_t
 * (the third bit is currently unused).
 */

/* Note that "0" is double precision,
 * so that you get that by default if
 * you forget a flag.
 */
#define GSL_PREC_DOUBLE  0
#define GSL_PREC_SINGLE  1
#define GSL_PREC_APPROX  2

#ifdef HAVE_INLINE
extern inline unsigned int
GSL_MODE_PREC(gsl_mode_t mt)
{ return  (mt & (unsigned int)7); }
#else
#define GSL_MODE_PREC(mt) ((mt) & (unsigned int)7)
#endif



/* Here are some predefined generic modes.
 */
#define GSL_MODE_DEFAULT  0


#endif /* !GSL_MODE_H */
