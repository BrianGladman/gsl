#ifndef GSL_SIMAN_H
#define GSL_SIMAN_H
#include <stdlib.h>

/* types for the function pointers passed to gsl_siman_solve */

typedef double (*gsl_Efunc_t) (void *xp);
typedef void (*gsl_siman_step_t) (const gsl_rng *r, void *xp, double step_size);
typedef double (*gsl_siman_metric_t) (void *xp, void *yp);
typedef void (*gsl_siman_print_t) (void *xp);

/* this structure contains all the information needed to structure the
   search, beyond the energy function, the step function and the
   initial guess. */

typedef struct {
  int n_tries;		/* how many points to try for each step */
  int iters_fixed_T;	/* how many iterations at each temperature? */
  double step_size;	/* max step size in the random walk */
  /* the following parameters are for the Boltzmann distribution */
  double k, t_initial, mu_t, t_min;
} gsl_siman_params_t;

struct CA_rule {
  char *str;
  int k, n, length;
};

typedef union u_Element {
  double D1;
  double D2[2];
  double D3[3];
  struct CA_rule rule;
} Element ;

/* prototype for the workhorse function */

void gsl_siman_Usolve(const gsl_rng * r, 
		      Element *x0_p, double (*Efunc)(Element x),
		      void (*take_step)(Element *x_p, double step_size),
		      double distance(Element x, Element y),
		      void print_position(Element x),
		      gsl_siman_params_t params);

void gsl_siman_solve(const gsl_rng * r, 
		     void *x0_p, gsl_Efunc_t Ef,
		     gsl_siman_step_t take_step,
		     gsl_siman_metric_t distance,
		     gsl_siman_print_t print_position,
		     size_t element_size,
		     gsl_siman_params_t params);

void 
gsl_siman_solve_many (const gsl_rng * r, void *x0_p, gsl_Efunc_t Ef,
		      gsl_siman_step_t take_step,
		      gsl_siman_metric_t distance,
		      gsl_siman_print_t print_position,
		      size_t element_size,
		      gsl_siman_params_t params);

#endif /* GSL_SIMAN_H */
