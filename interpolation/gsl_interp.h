/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_INTERP_H_
#define GSL_INTERP_H_


/* accelerator heuristics */
enum {
  GSL_INTERP_LOCALSTEP,
  GSL_INTERP_RANDOMACCESS
};


/* evaluation accelerator */
typedef struct {
  int cache_lo;
  int cache_hi;
  int cache_size;
  int heuristic;
  unsigned long miss_count;
  unsigned long hit_count;
}
gsl_interp_accel;


/* general interpolation object */
struct _gsl_interp_obj_struct {
  int     (*eval_impl)   (const struct _gsl_interp_obj_struct *, const double xa[], const double ya[], double x, gsl_interp_accel *, double * y);
  int     (*eval_na_impl) (const struct _gsl_interp_obj_struct *, const double xa[], const double ya[], double x, double * y);
  void    (*free)        (struct _gsl_interp_obj_struct *);
  double  xmin;
  double  xmax;
  int     size;
};
typedef  struct _gsl_interp_obj_struct  gsl_interp_obj;


/* interpolation object factory */
typedef struct {
  const char * name;
  gsl_interp_obj *  (*create) (const double x_array[], const double y_array[], int size);
}
gsl_interp_obj_factory;


/* available factories */
extern const gsl_interp_obj_factory   gsl_interp_linear_factory;
extern const gsl_interp_obj_factory   gsl_interp_cspline_natural_factory;
extern const gsl_interp_obj_factory   gsl_interp_cspline_periodic_factory;
extern const gsl_interp_obj_factory   gsl_interp_cspline_notnode_factory;
extern const gsl_interp_obj_factory   gsl_interp_cspline_fixed_factory;


gsl_interp_accel *
gsl_interp_accel_new(int heuristic, int cache_size);

void
gsl_interp_accel_free(gsl_interp_accel * a);


int
gsl_interp_eval_impl(const gsl_interp_obj * obj,
                     const double xa[], const double ya[], double x,
                     gsl_interp_accel * a, double * y
                     );

int
gsl_interp_eval_e(const gsl_interp_obj * obj,
                  const double xa[], const double ya[], double x,
                  gsl_interp_accel * a, double * y
                  );

double
gsl_interp_eval(const gsl_interp_obj * obj,
                const double xa[], const double ya[], double x,
                gsl_interp_accel * a
                );

int
gsl_interp_eval_noaccel_impl(const gsl_interp_obj * obj,
                             const double xa[], const double ya[], double x,
                             double * y
                             );

int
gsl_interp_eval_noaccel_e(const gsl_interp_obj * obj,
                          const double xa[], const double ya[], double x,
                          double * y
                          );

double
gsl_interp_eval_noaccel(const gsl_interp_obj * obj,
                        const double xa[], const double ya[], double x
                        );

void
gsl_interp_obj_free(gsl_interp_obj * interp_obj);


#endif  /* !GSL_INTERP_H_ */
