/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_POLY_H__
#define __GSL_SF_POLY_H__

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


/* c[0] + c[1] x + c[2] x^2 + ... + c[len-1] x^(len-1)
 *
 * exceptions: none
 */
double gsl_sf_poly_eval(const double c[], const int len, const double x);


#ifdef HAVE_INLINE
extern inline
double gsl_sf_poly_eval(const double c[], const int len, const double x)
{
  int i;
  double ans = c[len-1];
  for(i=len-1; i>0; i--) ans = c[i-1] + x * ans;
  return ans;
}
#endif /* HAVE_INLINE */


__END_DECLS

#endif /* __GSL_SF_POLY_H__ */
