/* Author: MJB */
/* RCS: $Id$ */
#include <gsl_rng.h>

enum {ESTIMATE_STYLE_NR = -1,  ESTIMATE_STYLE_CORRELATED_MC = 0,  
      ESTIMATE_STYLE_MC = 1};

extern unsigned long min_calls;
extern unsigned long min_calls_per_bisection;
extern double dither;
extern double estimate_frac;
extern double ALPHA;
extern int estimate_style;

typedef struct {
  unsigned long min_calls;
  unsigned long min_calls_per_bisection;
  double dither;
  double estimate_frac;
  double ALPHA;
  int estimate_style;
  int depth;
  int init_done;
  int check_done;
  gsl_rng *ranf;
} gsl_monte_miser_state; 

int gsl_monte_miser(const gsl_rng * r,
		    double (*func)(double []), double xl[], double xh[], 
		    unsigned long num_dim, unsigned long calls, 
		    double *ave, double *var);


gsl_monte_miser_state* gsl_monte_miser_alloc(void);

int gsl_monte_miser_validate(gsl_monte_miser_state* state,
			     double xl[], double xu[], 
			     unsigned long num_dim, unsigned long calls);

int gsl_monte_miser_init(gsl_monte_miser_state* state);


