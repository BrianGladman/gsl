/* Author: MJB */
/* RCS: $Id$ */

#ifndef __GSL_MONTE_MISER_H__
#define __GSL_MONTE_MISER_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>

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


__END_DECLS

#endif /* __GSL_MONTE_MISER_H__ */
