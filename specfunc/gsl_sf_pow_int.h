/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_POW_INT_H_
#define GSL_SF_POW_INT_H_

/* Calculate integer power. Especially good
   for small integer powers, for which it can
   be an order of magnitude faster than a
   call to standard library pow().
 */
double gsl_sf_pow_int(double x, int n);


#endif /* GSL_SF_POW_INT_H_ */
