/* Author:  B. Gough and G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_MACHINE_H
#define GSL_MACHINE_H

/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
#define GSL_MACH_EPS		1.0e-15
#define GSL_SQRT_MACH_EPS	3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       -34.54

/* machine limit constants */
#define GSL_SQRT_DBL_MIN        2.e-154
#define GSL_SQRT_DBL_MAX        5.e+153
#define GSL_ROOT3_DBL_MIN       3.42e-103
#define GSL_ROOT3_DBL_MAX       2.92e+102
#define GSL_ROOT4_DBL_MIN       1.414213e-77
#define GSL_ROOT4_DBL_MAX       7.071067e+76
#define GSL_LOG_DBL_MIN       	(DBL_MIN_10_EXP * M_LN10)
#define GSL_LOG_DBL_MAX       	(DBL_MAX_10_EXP * M_LN10)


#endif  /* !GSL_MACHINE_H */
