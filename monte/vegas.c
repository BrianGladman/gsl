/* Author: MJB */
/* RCS $Id$ */

/* This is an implementation of the adaptive Monte-Carlo algorithm
   of G. P. Lepage, originally described in J. Comp. Phys. 27, 192(1978).
   The current version of the algorithm was described in the Cornell
   preprint CLNS-80/447 of March, 1980.

   This code follows most closely the c version by D.R.Yennie, coded
   in 1984.

   The input coordinates are x[j], with upper and lower limits xu[j]
   and xl[j].  The integration length in the j-th direction is
   delx[j].  Each coordinate x[j] is mapped linearly to a variable
   y[j] in the range 0 to 1. This range is divided into bins with
   boundaries y_bin[i][j], where i=0 corresponds to y=0 and i=bins to
   y=1.  The integral contribution from the bin preceding y_bin[i][j]
   is called bin_sum[i][j] (thus bin_sum[0][j] is not used).  The grid
   is refined (ie, bins are adjusted) using grid_sum[i][j] which is
   some variation on the squared sum.  A third parameter used in
   defining the real coordinate using random numbers is called z.  It
   ranges from 0 to bins.  Its integer part gives the lower index of
   the bin into which a call is to be placed, and the remainder gives
   the location inside the bin.  

   The variable alpha controls how "stiff" the rebinning algorithm is.  
   alpha = 0 means never change the grid.  Alpha is typically set between
   1 and 2.

   The variable accuracy allows the program to terminate when the desired
   accuracty has been achieved.
*/

/* configuration headers */
#include <config.h>

/* standard headers */
#include <math.h>
#include <stdio.h>

/* gsl headers */
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_rng.h> 

/* lib-specific headers */

#define GSL_V_TINY 1.0e-30

#include <gsl_monte_vegas.h>
#include <vegas_print.h>

#define myMAX(a,b) ((a) >= (b) ? (a) : (b))

/* Setting the variable acc allows the user to terminate the computation
   when the desired accuracy has been achieved.
*/

/* All the state variables and workspace variables should be put into
   a structure or two, say vegas_flags and vegas_workspace.
*/

enum {MODE_IMPORTANCE = 1, MODE_IMPORTANCE_ONLY = 0, MODE_STRATIFIED = -1};


/* predeclare functions */

void adjust_bins(gsl_monte_vegas_state *state,
		 double bin[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM], 
		 double weight[GSL_V_BINS_MAX+1], 
		 double pts_per_bin, int j, int bins_prev, int bins);
inline int change_box_cord(int box_cord[GSL_V_MAX_DIM], 
			   int ng, unsigned long j);
inline void init_array(double array[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM], 
		       int imax, unsigned long jmax);


int gsl_monte_vegas_integrate(gsl_monte_vegas_state *state,
		    gsl_monte_f_T fxn, double xl[], double xu[], 
		    unsigned long num_dim, unsigned long calls,
		    double* tot_int, double* tot_sig, double* chi_sq)
{
  int i, j, k;
  int status;

  double cum_int = 1.0;
  double cum_sig = 1.0;

  gsl_rng *r;
  status = gsl_monte_vegas_validate(state, xl, xu, num_dim, calls);
  r = state->ranf;

  if (state->stage == 0) {
    init_array(state->y_bin, GSL_V_BINS_MAX, num_dim);
    for (j = 0; j < num_dim; ++j)
      state->y_bin[1][j] = 1.;
    state->bins_prev = 1;
    state->vol = 1;
    for (j = 0; j < num_dim; ++j) {
      state->delx[j] = xu[j] - xl[j];
      state->vol *= state->delx[j];
    }
    
    if (state->verbose >= 0 ) {
      vegas_open_log(state);
      print_lim(state, xl, xu, num_dim);
    }
  } /* stage == 0 */
  
  if (state->stage <= 1) {
    state->wtd_int_sum = 0;
    state->sum_wgts = 0;
    state->chi_sum = 0;
    state->it_num = 1;
  } /* stage <= 1 */

  if (state->stage <= 2) {
    state->bins = GSL_V_BINS_MAX;
    state->boxes = 1;
    if (state->mode != MODE_IMPORTANCE_ONLY) {
      state->boxes = floor( pow(calls/2.0, 1.0/num_dim) ); 
      state->mode = MODE_IMPORTANCE;
      if ((2 * state->boxes - GSL_V_BINS_MAX) >= 0) {
	state->mode = MODE_STRATIFIED;
	i = state->boxes / GSL_V_BINS_MAX;
	if ( state->boxes % GSL_V_BINS_MAX ) 
	  ++i;
	state->bins = state->boxes / i;
	state->boxes = i * state->bins;
      }
    }
    
    k = floor( pow( (double) state->boxes, (double) num_dim) );
    state->calls_per_box = myMAX(calls/k, 2);
    calls = state->calls_per_box * k;
    
    /* total volume of x-space/(avg num of calls/bin) */
    state->jac = state->vol *
      pow( (double) state->bins, (double) num_dim) / calls;
    
    /* If the number of bins changes from the previous invocation, bins
       are expanded or contracted accordingly, while preserving bin
       density */
    
    if (state->bins != state->bins_prev) {
      double weight[GSL_V_BINS_MAX+1], tot_weight;
      
      /* weight is ratio of bin sizes */
      tot_weight = (double) state->bins_prev / state->bins;
      for (i = 1; i <= state->bins_prev; ++i)
	weight[i] = 1;
      for (j = 0; j < num_dim; ++j)
	adjust_bins(state, state->y_bin, weight, tot_weight, j, 
		    state->bins_prev, state->bins);
      state->bins_prev = state->bins;
    }

    if (state->verbose >= 0) {
      print_head(state, 
		 num_dim, calls, state->it_num, state->bins, state->boxes);
    }
  } /* stage <= 2 */

  
  state->it_start = state->it_num;
  cum_int = 1.;
  cum_sig = 1.;
  for (; state->it_num <= state->max_it_num && 
	 state->acc*fabs(cum_int) < cum_sig; 
       ++state->it_num) {
    double intgrl, intgrl_sq, sig;
    double wgt, x[GSL_V_MAX_DIM];
    int    box_cord[GSL_V_MAX_DIM];
    for (j = 0; j < num_dim; ++j)
      box_cord[j] = 0;
    init_array(state->grid_sum, state->bins, num_dim);
    init_array(state->bin_sum, state->bins, num_dim);
    intgrl = 0.;
    sig = 0.;
    do {
      int    bin_cord[GSL_V_MAX_DIM];
      double f, f_sq, f_sum, f_sq_sum;

      for (k = 1, f_sum = 0., f_sq_sum = 0.; k <= state->calls_per_box; ++k) {
	double jacbin, y, z;

	jacbin = state->jac;
	for (j = 0; j < num_dim; ++j) {
	  double binwdth;

	  /* box_cord + ran gives the position in the box units, while
	     z is the position in bin units.
	  */
	  /*	  z = (box_cord[j] + gsl_rng_uniform(r) ) * bins_per_box; */
	  z = (box_cord[j] + gsl_rng_uniform(r) ) * state->bins / state->boxes; 
	  bin_cord[j] = z;
	  if (bin_cord[j] == 0) {
	    binwdth = state->y_bin[1][j];
	    y = z * binwdth;
	  } 
	  else {
	    binwdth = state->y_bin[bin_cord[j] + 1][j] - state->y_bin[bin_cord[j]][j];
	    y = state->y_bin[bin_cord[j]][j] + (z - bin_cord[j]) * binwdth;
	  }
	  x[j] = xl[j] + y * state->delx[j];
	  if (j > 2) {
	    /* fprintf(stderr, "z=%f,x[%d]=%f\n", z,j, x[j]);*/
	  }
	  jacbin *= binwdth;
	}
	/*
	 * At this point, jacbin is the real volume of the bin/(ave. no.
	 * calls/bin) 
	 */
	f = jacbin * fxn(x);
	
	f_sq = f * f;
	f_sum += f;
	f_sq_sum += f_sq;
	for (j = 0; j < num_dim; ++j) {
	  state->bin_sum[bin_cord[j] + 1][j] += f;
	  if (state->mode != MODE_STRATIFIED)
	    state->grid_sum[bin_cord[j] + 1][j] += f_sq;
	}
      } /* end of k loop */

      f_sq_sum = sqrt(f_sq_sum * state->calls_per_box);
      f_sq_sum = (f_sq_sum - f_sum) * (f_sq_sum + f_sum);
      intgrl += f_sum;
      if (f_sq_sum <= 0.0) f_sq_sum = GSL_V_TINY;
      sig += f_sq_sum;
      if (state->mode == MODE_STRATIFIED) {
	for (j = 0; j < num_dim; ++j)
	  state->grid_sum[bin_cord[j] + 1][j] += f_sq_sum;
      }
    } while ( change_box_cord ( box_cord, state->boxes, num_dim-1) );
    /* end of box_cord loop */

    /* Compute final results for this iteration   */

    sig = sig / (state->calls_per_box - 1);
    intgrl_sq = intgrl * intgrl;
    wgt = 1. / sig;
    state->wtd_int_sum += intgrl * wgt;
    state->sum_wgts += wgt;
    state->chi_sum += intgrl_sq * wgt;
    cum_int = state->wtd_int_sum / state->sum_wgts;
    *chi_sq = 0.;
    if (state->it_num != 1)
      *chi_sq = (state->chi_sum - state->wtd_int_sum * cum_int) / 
	(state->it_num - 1.);
    cum_sig = sqrt(1 / state->sum_wgts);

    if (state->verbose >= 0) {
      print_res(state, 
	      state->it_num, intgrl, sqrt(sig), cum_int, cum_sig, *chi_sq);
      if (state->it_num == state->max_it_num && state->verbose > 0)
	print_grid(state, num_dim);
    }

    /* Adjust the grid  */
    for (j = 0; j < num_dim; ++j) {
      double oldg, newg, grid_tot_j, tot_weight;
      double weight[GSL_V_BINS_MAX+1];

      oldg = state->grid_sum[1][j];
      newg = state->grid_sum[2][j];
      state->grid_sum[1][j] = (oldg + newg) / 2;
      grid_tot_j = state->grid_sum[1][j];

      /* This implements gs[i][j] = (gs[i-1][j]+gs[i][j]+gs[i+1][j])/3 */
      for (i = 2; i < state->bins; ++i) {
	state->grid_sum[i][j] = oldg + newg;
	oldg = newg;
	newg = state->grid_sum[i + 1][j];
	state->grid_sum[i][j] = (state->grid_sum[i][j] + newg) / 3;
	grid_tot_j += state->grid_sum[i][j];
      }
      state->grid_sum[state->bins][j] = (newg + oldg) / 2;

      grid_tot_j += state->grid_sum[state->bins][j];
      for (i = 1, tot_weight = 0; i <= state->bins; ++i) {
	weight[i] = 0;
	if (state->grid_sum[i][j] > 0) {
	  oldg = grid_tot_j / state->grid_sum[i][j];
	  /* damped change */
	  weight[i] = pow(((oldg - 1) / oldg / log(oldg)), state->alpha);
	}
	tot_weight += weight[i];
      }
      adjust_bins(state, state->y_bin, weight, tot_weight/state->bins, j, 
		  state->bins, state->bins);
    }
  } /* for it_num */

  *tot_int = cum_int;
  *tot_sig = cum_sig;

  if ( (state->it_num > state->max_it_num) || 
       (state->acc*fabs(cum_int) > cum_sig) ) {
    /* should throw an error of some kind */
  }

  /* final stuff */
  if (state->verbose >= 0 ) {
    vegas_close_log(state);
  }

  return status;
}

/* Adjust bin boundaries so that each bin has the same number of
   points/bin.  Do this by stepping through the current bins, moving
   the boundary as we go and counting the number of points until we
   exceed the desired number.  Then we move the boundary back an
   appropriate amount.  
*/

void adjust_bins(gsl_monte_vegas_state *state,
		 double bin[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM], 
		 double weight[GSL_V_BINS_MAX+1], 
		 double pts_per_bin, int j, int bins_prev, int bins)
{
  int    i, k;
  double xold, xnew, xin[GSL_V_BINS_MAX+1], dw;

  xnew = 0;
  dw = 0;
  i = 1;
  for (k = 1; k <= bins_prev; ++k) {
    dw += weight[k];
    xold = xnew;
    xnew = bin[k][j];
    for (; dw > pts_per_bin; ++i) {
      dw -= pts_per_bin;
      xin[i] = xnew - (xnew - xold) * dw / weight[k];
    }
  }
  for (i = 1; i < bins; ++i)
    bin[i][j] = xin[i];
  bin[bins][j] = 1;
  return;
}

inline int change_box_cord(int box_cord[GSL_V_MAX_DIM], 
			   int ng, unsigned long j_start)
{
  int j = j_start;
  while ( j >= 0 ) {
    ++box_cord[j];
    if ( (box_cord[j] %= ng) != 0) 
      return (1);
    j--;
  }
  return (0);
}

inline void init_array(double array[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM], 
		       int imax, unsigned long jmax)
{
  int i, j;

  for (j = 0; j < jmax; ++j) {
    for (i = 0; i <= imax; ++i)
      array[i][j] = 0;
  }
}

gsl_monte_vegas_state* gsl_monte_vegas_alloc(void)
{
  gsl_monte_vegas_state *s =  
    (gsl_monte_vegas_state *) malloc(sizeof (gsl_monte_vegas_state));
  
  if ( s == (gsl_monte_vegas_state*) NULL) {
    GSL_ERROR_RETURN ("failed to allocate space for miser state struct",
                        GSL_ENOMEM, 0);
  }

  return s;
}

int gsl_monte_vegas_validate(gsl_monte_vegas_state* state,
			     double xl[], double xu[], 
			     unsigned long num_dim, unsigned long calls)
{
  unsigned long i;
  char warning[100];

  if (state == (gsl_monte_vegas_state*) NULL) {
    GSL_ERROR("Allocate state structure before calling!", GSL_EINVAL);

  }

  if (state->check_done) 
    return GSL_SUCCESS;
    
  if (num_dim <= 0) {
    sprintf(warning, "number of dimensions must be greater than zero, not %lu",
	    num_dim);
    GSL_ERROR(warning, GSL_EINVAL);
  }
  
  for (i=0; i < num_dim; i++ ) {
    if (xu[i] - xl[i] <= 0 ) {
      sprintf(warning, "xu[%lu] must be greater than xu[%lu]", i, i);
    GSL_ERROR(warning, GSL_EINVAL);
    }
    if (xu[i] - xl[i] > DBL_MAX) {
      sprintf(warning, 
	      "Range of integration is too large for cord %lu, please rescale", 
	      i);
      GSL_ERROR(warning, GSL_EINVAL);
    }
  }

  if ( calls <= 0 ) {
    sprintf(warning, "number of calls must be greater than zero, not %lu",
	    calls);
    GSL_ERROR(warning, GSL_EINVAL);
  }
  
  state->check_done = 1;

  return GSL_SUCCESS;
}  

/* Set some default values and whatever */
int gsl_monte_vegas_init(gsl_monte_vegas_state* state)
{

  if (state == (gsl_monte_vegas_state*) NULL) {
    GSL_ERROR("Allocate state structure before calling!", GSL_EFAULT);
  }

  state->stage = 0;
  state->acc = -1;
  state->alpha = 1.5;
  state->verbose = -1;
  state->max_it_num = 5;
  state->ranf = gsl_rng_alloc(gsl_rng_env_setup());

  state->init_done = 1;
  return GSL_SUCCESS;
}

void gsl_monte_vegas_free (gsl_monte_vegas_state* s)
{
  if (s == (gsl_monte_vegas_state*) NULL )
    GSL_ERROR_RETURN_NOTHING("Attempt to free null pointer", GSL_EFAULT);

  free (s);
}
