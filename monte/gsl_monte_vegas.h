/* header for the gsl "vegas" routines.  Mike Booth, May 1998 */
/* RCS $Id$ */

#ifndef GSL_MONTE_VEGAS_H
#define GSL_MONTE_VEGAS_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <stdio.h>

/* This will go away soon. */
#define GSL_V_BINS_MAX 50  /* even integer because will be divided by two. */
#define GSL_V_MAX_DIM 10

enum {GSL_VEGAS_MODE_IMPORTANCE = 1, 
      GSL_VEGAS_MODE_IMPORTANCE_ONLY = 0, 
      GSL_VEGAS_MODE_STRATIFIED = -1};

typedef struct {
  /* control variables */
  double acc;
  double alpha;
  int mode;
  int verbose;
  int max_it_num;
  int stage;

  /* state variables */
  int it_start;
  int bins_prev;
  int calls_per_box;
  int it_num;
  unsigned long num_dim;
  int bins;
  int boxes; /* boxes and bins are counted along the axes */
  int init_done;
  int check_done;
  gsl_rng* ranf;
  FILE* ostream;

  /* scratch variables preserved between calls to vegas1/2/3  */
  double jac;
  double wtd_int_sum; 
  double sum_wgts;
  double chi_sum;
  double vol;

  /* workspace */
  double delx[GSL_V_MAX_DIM];
  double grid_sum[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM];
  double bin_sum[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM];
  double y_bin[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM];

} gsl_monte_vegas_state;


int gsl_monte_vegas_integrate(gsl_monte_vegas_state *state,
		    gsl_monte_f_T fxn, double xl[], double xu[], 
		    unsigned long num_dim, unsigned long calls,
		    double* tot_int, double* tot_sig, double* chi_sq);

gsl_monte_vegas_state* gsl_monte_vegas_alloc(size_t num_dim);

int gsl_monte_vegas_validate(gsl_monte_vegas_state* state,
			     double xl[], double xu[], 
			     unsigned long num_dim, unsigned long calls);

int gsl_monte_vegas_init(gsl_monte_vegas_state* state);

void gsl_monte_vegas_free (gsl_monte_vegas_state* s);

#endif /* !GSL_MONTE_VEGAS_H */

