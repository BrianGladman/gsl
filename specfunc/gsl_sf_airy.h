/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_AIRY_H_
#define GSL_AIRY_H_

/* Airy function Ai(x) 
 * No error condition can occur.
 */
double gsl_sf_airy_Ai(double x);

/* Airy function Bi(x) with POSIX-ish error status return.
 * Returns GSL_EOVRFLW on overflow, GSL_SUCCESS on success.
 */
int gsl_sf_airy_Bi_e(double x, double *result);

/* Airy function Bi(x). Overflow can occur.
 */
double gsl_sf_airy_Bi(double x);

/* Airy function Bi(x) with exponential
 * prefactor removed when x > 0
 * No error condition can occur.
 */
double gsl_sf_airy_Bi_scaled(double x);


#endif /* GSL_AIRY_H_ */
