/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include "gsl_sf_poly.h"

/* checked OK [GJ] Tue May  5 12:19:56 MDT 1998 */
inline double gsl_sf_poly_eval(const double c[], const int len, const double x)
{
  int i;
  double ans = c[len-1];
  for(i=len-1; i>0; i--) ans = c[i-1] + x * ans;
  return ans;
}
