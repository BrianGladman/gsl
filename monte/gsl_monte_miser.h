/* Author: MJB */
/* RCS: $Id$ */

#ifndef GSL_MONTE_MISER_H
#define GSL_MONTE_MISER_H

#include <gsl_rng.h>
#include <gsl_monte.h>
#include <gsl_monte_plain.h>

enum {ESTIMATE_STYLE_NR = -1,  ESTIMATE_STYLE_CORRELATED_MC = 0,  
      ESTIMATE_STYLE_MC = 1};

typedef struct {
  unsigned long min_calls;
  unsigned long min_calls_per_bisection;
  double dither;
  double estimate_frac;
  double alpha;
  size_t num_dim;
  int estimate_style;
  int depth;
  int verbose;
  int init_done;
  int check_done;
  FILE* ostream;
  gsl_rng *ranf;
  gsl_monte_plain_state* plain_state;
} gsl_monte_miser_state; 

int gsl_monte_miser_integrate(gsl_monte_miser_state* state,
			      gsl_monte_f_T func, double xl[], double xh[], 
			      unsigned long num_dim, unsigned long calls, 
			      double *ave, double *var);


gsl_monte_miser_state* gsl_monte_miser_alloc(size_t num_dim);

int gsl_monte_miser_validate(gsl_monte_miser_state* state,
			     double xl[], double xu[], 
			     unsigned long num_dim, unsigned long calls);

int gsl_monte_miser_init(gsl_monte_miser_state* state);

void gsl_monte_miser_free(gsl_monte_miser_state* state);


#endif /* GSL_MONTE_MISER_H */
