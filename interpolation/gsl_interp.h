/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_INTERP_H_
#define GSL_INTERP_H_


/* iterator heuristics */
enum {
  GSL_INTERP_LOCALSTEP,
  GSL_INTERP_RANDOMACCESS
};

/* spline types */
enum {
  GSL_INTERP_CSPLINE_FREE,
  GSL_INTERP_CSPLINE_FIXED
};

/* iterator */
typedef struct {
  int cache_lo;
  int cache_hi;
  int cache_size;
  int type;
  unsigned long miss_count;
  unsigned long hit_count;
}
gsl_interp_iter;

/* linear interpolation data */
typedef struct {
  double  xmin;
  double  xmax;
  int     size;
}
gsl_interp_linear;

/* cubic spline interpolation data */
typedef struct {
  double    xmin;
  double    xmax;
  int       size;
  int       type;
  double *  y2;
}
gsl_interp_cspline;


gsl_interp_iter *
gsl_interp_iter_new(int type, int cache_size);

void
gsl_interp_iter_free(gsl_interp_iter * iter);


gsl_interp_linear *
gsl_interp_linear_new(const double x_array[], const double y_array[], int size);

int
gsl_interp_linear_eval_impl(const gsl_interp_linear * interp_lin,
                            gsl_interp_iter * iter,
                            const double x_array[], const double y_array[],
		            double x,
			    double * y
		            );

int
gsl_interp_linear_eval_direct_impl(const gsl_interp_linear * interp_lin,
                                   const double x_array[], const double y_array[],
		                   double x,
			           double * y
		                   );

void
gsl_interp_linear_free(gsl_interp_cspline * interp_cs);


gsl_interp_cspline *
gsl_interp_cspline_new(const double x[], const double y[], int size, int type);

int
gsl_interp_cspline_eval_impl(const gsl_interp_cspline * interp_cs,
                             gsl_interp_iter * iter,
                             const double x_array[], const double y_array[],
		             double x,
			     double * y
		             );

void
gsl_interp_cspline_free(gsl_interp_cspline * interp_cs);


#endif  /* !GSL_INTERP_H_ */
