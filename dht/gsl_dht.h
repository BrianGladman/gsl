/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_DHT_H
#define GSL_DHT_H


struct gsl_dht_transform_struct {
  size_t    size;  /* size of the sample arrays to be transformed    */
  double    nu;    /* Bessel function order                          */
  double    xmax;  /* the upper limit to the x-sampling domain       */
  double    kmax;  /* the upper limit to the k-sampling domain       */
  double *  j;     /* array of computed J_nu zeros, j_{nu,s} = j[s]  */
  double *  Jjj;   /* transform numerator, J_nu(j_i j_m / j_N)       */
  double *  J2;    /* transform denominator, J_{nu+1}^2(j_m)         */
};
typedef struct gsl_dht_transform_struct gsl_dht_transform;


/* Create a new transform object for a given size
 * sampling array on the domain [0, xmax].
 */
gsl_dht_transform * gsl_dht_transform_new(size_t size, double nu, double xmax);


/* Recalculate a transform object for given values of nu, xmax.
 * You cannot change the size of the object since the internal
 * allocation is reused.
 */
int gsl_dht_transform_recalc_impl(gsl_dht_transform * t, double nu, double xmax);
int gsl_dht_transform_recalc_e(gsl_dht_transform * t, double nu, double xmax);


/* The n'th computed x sample point for a given transform. */
double gsl_dht_transform_x_sample(const gsl_dht_transform * t, int n);


/* The n'th computed k sample point for a given transform. */
double gsl_dht_transform_k_sample(const gsl_dht_transform * t, int n);


/* Free a transform object. */
void gsl_dht_transform_free(gsl_dht_transform * t);


/* Perform a transform on a sampled array. */
int gsl_dht_transform_apply(const gsl_dht_transform * t, double * f_in, double * f_out);


#endif  /* !GSL_DHT_H */
