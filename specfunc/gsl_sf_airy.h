/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_AIRY_H_
#define GSL_AIRY_H_

/* Airy function Ai(x) */
double gsl_sf_airy_Ai(double x);

/* Airy function Bi(x) */
double gsl_sf_airy_Bi(double x);

/* Airy function Bi(x) with exponential
   prefactor removed when x > 0
 */
double gsl_sf_airy_Bi_scaled(double x);


#endif /* GSL_AIRY_H_ */
