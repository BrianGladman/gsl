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
#include <gsl_ran.h>

/* lib-specific headers */

#define GSL_V_TINY 1.0e-30

#include "gsl_vegas.h"
#include "gsl_vegas_print.h"

#define myMAX(a,b) ((a) >= (b) ? (a) : (b))

/* Setting the variable acc allows the user to terminate the computation
   when the desired accuracy has been achieved.
*/

/* All the state variables and workspace variables should be put into
   a structure or two, say vegas_flags and vegas_workspace.
*/

enum {MODE_IMPORTANCE = 1, MODE_IMPORTANCE_ONLY = 0, MODE_STRATIFIED = -1};

int    calls = 1000;
double acc = -1, alpha = 1.5;
int    mode, verbose = -1, max_it_num = 5;

static int    it_start, bins_prev, calls_per_box, it_num, bins, boxes;
static double delx[GSL_V_MAX_DIM], grid_sum[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM];
static double bin_sum[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM];
static double y_bin[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM];

static double jac, wtd_int_sum, sum_wgts, chi_sum, vol;

/* predeclare functions */

void adjust_bins(double bin[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM], 
		 double weight[GSL_V_BINS_MAX+1], 
		 double pts_per_bin, int j, int n1, int n2);
inline int change_box_cord(int box_cord[GSL_V_MAX_DIM], int ng, int j);
inline void init_array(double array[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM], 
		       int n1, int n2);


int gsl_monte_vegas(gsl_monte_f_T fxn, double xl[], double xu[], int num_dim,
		    double* tot_int, double* tot_sig, double* chi_sq_ptr)
{
  int j;
  int status;

  init_array(y_bin, GSL_V_BINS_MAX, num_dim);
  for (j = 0; j < num_dim; ++j)
    y_bin[1][j] = 1.;
  bins_prev = 1;
  vol = 1;
  for (j = 0; j < num_dim; ++j) {
    delx[j] = xu[j] - xl[j];
    vol *= delx[j];
  }

  if (verbose >= 0 ) {
    vegas_open_log();
    prn_lim(xl, xu, num_dim);
  }
  status = gsl_monte_vegas1(fxn, xl, xu, num_dim, tot_int, tot_sig, chi_sq_ptr);

  if (verbose >= 0 ) {
    vegas_close_log();
  }

  return status;
}

int gsl_monte_vegas1(gsl_monte_f_T fxn, double xl[], double xu[], int num_dim,
		     double* tot_int, double* tot_sig, double* chi_sq_ptr)
{
  int status;

  wtd_int_sum = 0;
  sum_wgts = 0;
  chi_sum = 0;
  it_num = 1;

  status = gsl_monte_vegas2(fxn, xl, xu, num_dim, tot_int, tot_sig, chi_sq_ptr);
  return status;
}

int gsl_monte_vegas2(gsl_monte_f_T fxn, double xl[], double xu[], int num_dim,
		     double* tot_int, double* tot_sig, double* chi_sq_ptr)
{

  /* At this point, the calculation is divided into three possible
     cases: mode = MODE_IMPORTANCE_ONLY: obtained only by user choice
     at the beginning.  There is only one box, and the number of bins
     is set to GSL_V_BINS_MAX.  The points are distributed over the whole
     grid randomly, but the bins adjust to emphasize the regions where
     the function is largest.  mode = MODE_IMPORTANCE or
     MODE_STRATIFIED is decided by a test - if it is not possible to
     have at least 2 bins/box with the maximal number of boxes, then
     MODE_STRATIFIED is chosen.  Generally, this happens for low
     dimensions where the vegas algorithm has little to offer.
     Choosing MODE_IMPORTANCE concentrates increments where the
     integrand is largest and MODE_STATIFIED where the error is
     largest.

     In addition to the bin structure, the volume is divided into
     boxes, with a equal number of calls (calls_per_box) going into
     each box.  In both cases, the number of calls is given by calls =
     (calls_per_box * boxes)^num_dim and because we round things to
     integers is usually a bit smaller than the nominal number
     specified.  For MODE_IMPORTANCE, the number of bins in each
     dimension is GSL_V_BINS_MAX, while the number of boxes in each
     dimension is smaller (because of the above-mentioned
     restriction).  Thus, each box contains several (not necessarily
     an integral number) bins.  For MODE_STRATIFIED, bins is an
     integer multiple of boxes, and bins is smaller than GSL_V_BINS_MAX.  */

  int    i, j, k;
  int status;

  bins = GSL_V_BINS_MAX;
  boxes = 1;
  if (mode != MODE_IMPORTANCE_ONLY) {
    boxes = floor( pow(calls/2.0, 1.0/num_dim) ); 
    mode = MODE_IMPORTANCE;
    if ((2 * boxes - GSL_V_BINS_MAX) >= 0) {
      mode = MODE_STRATIFIED;
      i = boxes / GSL_V_BINS_MAX + 1;
      bins = boxes / i;
      boxes = i * bins;
    }
  }

  k = floor( pow( (double) boxes, (double) num_dim) );
  calls_per_box = myMAX(calls/k, 2);
  calls = calls_per_box * k;
  /*  fprintf(stderr, "mode=%d, calls_per_box=%d, bins=%d, boxes=%d\n", 
	  mode, calls_per_box, bins, boxes);
  */
  /* bins_per_box = (double) bins / boxes; */

  /* total volume of x-space/(avg num of calls/bin) */
  jac = vol * pow( (double) bins, (double) num_dim) / calls;

  /* If the number of bins changes from the previous invocation, bins
     are expanded or contracted accordingly, while preserving bin
     density */

  if (bins != bins_prev) {
    double weight[GSL_V_BINS_MAX+1], tot_weight;

    tot_weight = (double) bins_prev / bins;	/* ratio of bin sizes */
    for (i = 1; i <= bins_prev; ++i)
      weight[i] = 1;
    for (j = 0; j < num_dim; ++j)
      adjust_bins(y_bin, weight, tot_weight, j, bins_prev, bins);
    bins_prev = bins;
  }
  if (verbose >= 0) {
    prn_head(num_dim, calls, it_num, max_it_num, acc, verbose, alpha, mode, 
	     bins, boxes);
  }
  status = gsl_monte_vegas3(fxn, xl, xu, num_dim, tot_int, tot_sig, chi_sq_ptr);
  return status;
}

int gsl_monte_vegas3(gsl_monte_f_T fxn, double xl[], double xu[], int num_dim,
		     double* tot_int, double* tot_sig, double* chi_sq_ptr)
{

  /* Start of iteration  */

  int    i, j, k;
  int status = 0;

  double cum_int, cum_sig;
  double  chi_sq = 0;

  it_start = it_num;
  cum_int = 1.;
  cum_sig = 1.;
  for (; it_num <= max_it_num && acc*fabs(cum_int) < cum_sig; ++it_num) {
    double intgrl, intgrl_sq, sig;
    double wgt, x[GSL_V_MAX_DIM];
    int    box_cord[GSL_V_MAX_DIM];
    for (j = 0; j < num_dim; ++j)
      box_cord[j] = 0;
    init_array(grid_sum, bins, num_dim);
    init_array(bin_sum, bins, num_dim);
    intgrl = 0.;
    sig = 0.;
    do {
      int    bin_cord[GSL_V_MAX_DIM];
      double f, f_sq, f_sum, f_sq_sum;

      for (k = 1, f_sum = 0., f_sq_sum = 0.; k <= calls_per_box; ++k) {
	double jacbin, y, z;

	jacbin = jac;
	for (j = 0; j < num_dim; ++j) {
	  double binwdth;

	  /* box_cord + ran gives the position in the box units, while
	     z is the position in bin units.
	  */
	  /*	  z = (box_cord[j] + gsl_ran_uniform() ) * bins_per_box; */
	  z = (box_cord[j] + gsl_ran_uniform() ) * bins / boxes; 
	  bin_cord[j] = z;
	  if (bin_cord[j] == 0) {
	    binwdth = y_bin[1][j];
	    y = z * binwdth;
	  } 
	  else {
	    binwdth = y_bin[bin_cord[j] + 1][j] - y_bin[bin_cord[j]][j];
	    y = y_bin[bin_cord[j]][j] + (z - bin_cord[j]) * binwdth;
	  }
	  x[j] = xl[j] + y * delx[j];
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
	  bin_sum[bin_cord[j] + 1][j] += f;
	  if (mode != MODE_STRATIFIED)
	    grid_sum[bin_cord[j] + 1][j] += f_sq;
	}
      } /* end of k loop */

      f_sq_sum = sqrt(f_sq_sum * calls_per_box);
      f_sq_sum = (f_sq_sum - f_sum) * (f_sq_sum + f_sum);
      intgrl += f_sum;
      if (f_sq_sum <= 0.0) f_sq_sum = GSL_V_TINY;
      sig += f_sq_sum;
      if (mode == MODE_STRATIFIED) {
	for (j = 0; j < num_dim; ++j)
	  grid_sum[bin_cord[j] + 1][j] += f_sq_sum;
      }
    } while ( change_box_cord ( box_cord, boxes, num_dim-1) );
    /* end of box_cord loop */

    /* Compute final results for this iteration   */

    sig = sig / (calls_per_box - 1);
    intgrl_sq = intgrl * intgrl;
    wgt = 1. / sig;
    wtd_int_sum += intgrl * wgt;
    sum_wgts += wgt;
    chi_sum += intgrl_sq * wgt;
    cum_int = wtd_int_sum / sum_wgts;
    chi_sq = 0.;
    if (it_num != 1)
      chi_sq = (chi_sum - wtd_int_sum * cum_int) / (it_num - 1.);
    cum_sig = sqrt(1 / sum_wgts);

    if (verbose >= 0) {
      prn_res(it_num, intgrl, sqrt(sig), cum_int, cum_sig, chi_sq);
      if (it_num == max_it_num && verbose > 0)
	prn_grid(y_bin, bin_sum, num_dim, bins, verbose);
    }

    /* Adjust the grid  */
    for (j = 0; j < num_dim; ++j) {
      double oldg, newg, grid_tot_j, weight[GSL_V_BINS_MAX+1], tot_weight;

      oldg = grid_sum[1][j];
      newg = grid_sum[2][j];
      grid_sum[1][j] = (oldg + newg) / 2;
      grid_tot_j = grid_sum[1][j];

      /* This implements gs[i][j] = (gs[i-1][j]+gs[i][j]+gs[i+1][j])/3 */
      for (i = 2; i < bins; ++i) {
	grid_sum[i][j] = oldg + newg;
	oldg = newg;
	newg = grid_sum[i + 1][j];
	grid_sum[i][j] = (grid_sum[i][j] + newg) / 3;
	grid_tot_j += grid_sum[i][j];
      }
      grid_sum[bins][j] = (newg + oldg) / 2;

      grid_tot_j += grid_sum[bins][j];
      for (i = 1, tot_weight = 0; i <= bins; ++i) {
	weight[i] = 0;
	if (grid_sum[i][j] > 0) {
	  oldg = grid_tot_j / grid_sum[i][j];
	  /* damped change */
	  weight[i] = pow(((oldg - 1) / oldg / log(oldg)), alpha);
	}
	tot_weight += weight[i];
      }
      adjust_bins(y_bin, weight, tot_weight/bins, j, bins, bins);
    }
  } /* for it_num */

  *tot_int = cum_int;
  *tot_sig = cum_sig;
  *chi_sq_ptr = chi_sq;

  if ( (it_num > max_it_num) || (acc*fabs(cum_int) > cum_sig) ) {
    /* should throw an error of some kind */
  }
  return status;

}

/* Adjust bin boundaries so that each bin has the same number of
   points/bin.  Do this by stepping through the current bins, moving
   the boundary as we go and counting the number of points until we
   exceed the desired number.  Then we move the boundary back an
   appropriate amount.  
*/

void adjust_bins(double bin[GSL_V_BINS_MAX+1][GSL_V_MAX_DIM], 
		 double weight[GSL_V_BINS_MAX+1], 
		 double pts_per_bin, int j, int n1, int n2)
{
  int    i, k;
  double xold, xnew, xin[GSL_V_BINS_MAX+1], dw;

  xnew = 0;
  dw = 0;
  i = 1;
  for (k = 1; k <= n1; ++k) {
    dw += weight[k];
    xold = xnew;
    xnew = bin[k][j];
    for (; dw > pts_per_bin; ++i) {
      dw -= pts_per_bin;
      xin[i] = xnew - (xnew - xold) * dw / weight[k];
    }
  }
  for (i = 1; i < n2; ++i)
    bin[i][j] = xin[i];
  bin[n2][j] = 1;
  return;
}

inline int change_box_cord(int box_cord[GSL_V_MAX_DIM], int ng, int j_start)
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
		       int n1, int n2)
{
  int i, j;

  for (j = 0; j < n2; ++j) {
    for (i = 0; i <= n1; ++i)
      array[i][j] = 0;
  }
}
