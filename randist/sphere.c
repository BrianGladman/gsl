#include <config.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

#define USE_KNOP_SPHERE_METHOD 1 /* it really is better! */
#define USE_TRIG_FOR_CIRCLE    0 /* it's faster on Pentium, not on Sun */
#define USE_VON_NEUMANN_TRICK  1 /* avoids a sqrt() */

void
gsl_ran_dir_2d (const gsl_rng * r, double *x, double *y)
{
#if USE_TRIG_FOR_CIRCLE
    /* This is the obvious solution... */
    /* It ain't clever, but since sin/cos are often hardware accelerated,
     * it can be faster -- it is on my home Pentium -- than von Neumann's
     * solution below, or slower -- as it is on my Sun Sparc 20 at work
     */     
    double t;
    t = 6.2831853071795864*gsl_rng_uniform(r);  /* 2*PI */
    *px = cos(t);
    *py = sin(t);
#else /* dont USE_TRIG_FOR_CIRCLE */
    /* Avoids trig, but it does take an average of 8/pi = 2.55 calls to
     * the RNG, instead of one, as above.
     */
    double u,v,s;
    do {
        u = -1 + 2*gsl_rng_uniform(r);
        v = -1 + 2*gsl_rng_uniform(r);
        s = u*u+v*v;
    } while (s > 1.0 || s == 0.0);
#if USE_VON_NEUMANN_TRICK 
    /* See Knuth, v2, 3rd ed, p140 (exercise 23).
     * Note, no sin, cos, or sqrt !
     */
    *x = (u*u-v*v)/s;
    *y = 2*u*v/s;
#else  /* don't USE_VON_NEUMANN_TRICK */  
    /* Here is the more straightforward approach */
    /* Fewer total operations, but one of them is a sqrt */
    s = sqrt(s);
    *x = u/s;
    *y = v/s;
#endif  /* end USE_VON_NEUMANN_TRICK */
#endif  /* end USE_TRIG_FOR_CIRCLE */
}

void
gsl_ran_dir_3d (const gsl_rng * r, double *x, double *y, double *z)
{
    double s,a;
#if USE_KNOP_SPHERE_METHOD
    /* This is a variant of the algorithm for computing a random point
     * on the unit sphere; the algorithm is suggested in Knuth, v2,
     * 3rd ed, p136; and attributed to Robert E Knop, CACM, 13 (1970),
     * 326.
     */
    
    /* Begin with the polar method for getting x,y inside a unit circle
     */
    do {
        *x = -1 + 2*gsl_rng_uniform(r);
        *y = -1 + 2*gsl_rng_uniform(r);
        s = (*x)*(*x)+(*y)*(*y);
    } while (s > 1.0 || s == 0.0);

    *z = -1 + 2*s;              /* z uniformly distributed from -1 to 1 */
    a = 2*sqrt(1-s);            /* factor to adjust x,y so that x^2+y^2
                                 * is equal to 1-z^2 */
    *x *= a;
    *y *= a;
#else /* don't USE_KNOP_SPHERE_METHOD */
    /* Here is the more straightforward method: expect 18/pi=5.7
     * calls to RNG on average, instead of 2.55 with Knop method.
     */
    do {
        *x = -1 + 2*gsl_rng_uniform(r);
        *y = -1 + 2*gsl_rng_uniform(r);
        *z = -1 + 2*gsl_rng_uniform(r);
        s = (*x)*(*x)+(*y)*(*y)+(*z)*(*z);
    } while (s > 1.0 || s == 0.0);
    s = sqrt(s);
    *x /= s;
    *y /= s;
    *z /= s;
#endif  /* end USE_KNOP_SPHERE_METHOD */
}

void
gsl_ran_dir_nd (const gsl_rng * r, int n, double *x)
{
  double d;
  int i;
  /* See Knuth, v2, 3rd ed, p135-136.  The method is attributed to
   * G. W. Brown, in Modern Mathematics for the Engineer (1956).
   * The idea is that gaussians G(x) have the property that
   * G(x)G(y)G(z)G(...) is radially symmetric, a function only
   * r = sqrt(x^2+y^2+...)
   */
  d = 0;
  do {
      for (i=0; i<n; ++i) {
          x[i] = gsl_ran_gaussian(r,1.0);
          d += x[i]*x[i];
      }
  }
  while (d == 0);
  d = sqrt(d);
  for (i=0; i<n; ++i) {
      x[i] /= d;
  }
}      
