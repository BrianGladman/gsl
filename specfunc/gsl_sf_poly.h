/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_POLY_H_
#define GSL_SF_POLY_H_

/* c[0] + c[1] x + c[2] x^2 + ... + c[len-1] x^(len-1)
 */
double gsl_sf_poly_eval(const double c[], const int len, const double x);

#ifdef HAVE_INLINE
extern inline double gsl_sf_poly_eval(const double c[], const int len, const double x)
{
  int i;
  double ans = c[len-1];
  for(i=len-1; i>0; i--) ans = c[i-1] + x * ans;
  return ans;
}
#endif /* HAVE_INLINE */

#endif  /* !GSL_SF_POLY_H_ */
