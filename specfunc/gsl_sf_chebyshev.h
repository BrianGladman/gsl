#ifndef GSL_SF_CHEBYSHEV_H_
#define GSL_SF_CHEBYSHEV_H_


/* data for a Chebyshev series over a given interval */
struct gsl_sf_ChebSeries {
  double * c;   /* coefficients */
  int order;    /* order of expansion */
  double a;     /* lower interval point */
  double b;     /* upper interval point */
};


/* calculate a Chebyshev series of specified order over
   a specified interval, for a given function
 */
extern struct gsl_sf_ChebSeries * cheb_new(double (*func)(double),
					   double a, double b,
					   int order);

/* evaluate a Chebyshev series at a given point */
extern double gsl_sf_cheb_eval(double x, const struct gsl_sf_ChebSeries * cs);


/* free a Chebyshev series previously calculated with gsl_sf_cheb_new() */
extern void gsl_sf_cheb_free(struct gsl_sf_ChebSeries * cs);


#endif /* GSL_SF_CHEBYSHEV_H_ */
