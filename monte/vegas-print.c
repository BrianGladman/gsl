/*   Set of programs for printing vegas results   */
/* RCS $Id$ */

#include <config.h>
#include <stdio.h>
#include <gsl_monte_vegas.h>
#include <vegas_print.h>
/* predeclare functions */

/* The following variables are needed for the print functions. */

int vegas_open_log(gsl_monte_vegas_state* state)
{
  state->ostream = fopen("vegas.out", "a+");
  if (state->ostream == (FILE*) NULL)
    /* FIXME: Maybe better error code? */
    GSL_ERROR("failed to open stream on vegas.out", GSL_EFAILED);
  return 0;
}

int vegas_close_log(gsl_monte_vegas_state* state)
{
  if (state->ostream == (FILE*) NULL) 
    GSL_ERROR("attempted to close null file pointer", GSL_EFAULT);

  return fclose(state->ostream);
  
}

void print_lim(gsl_monte_vegas_state* state, 
	       double xl[], double xu[], unsigned long dim)
{
  int j;

  if (state->ostream == (FILE*) NULL) 
    return;

  fprintf(state->ostream, "The limits of integration are:\n");
  for (j = 0; j < dim; ++j)
    fprintf(state->ostream, "\nxl[%d]=%f    xu[%d]=%f", j, xl[j], j, xu[j]);
  fprintf(state->ostream, "\n");
  fflush(state->ostream);
}

void print_head(gsl_monte_vegas_state* state, 
		unsigned long num_dim, unsigned long calls, 
		int it_num, int bins, int boxes)
{
  if (state->ostream == (FILE*) NULL) 
    return;

  fprintf(state->ostream, 
	  "\nnum_dim=%lu, calls=%lu, it_num=%d, max_it_num=%d, acc=%.3f, ",
	  num_dim, calls, it_num, state->max_it_num, state->acc);
  fprintf(state->ostream, "verb=%d, alph=%.2f,\nmode=%d, bins=%d, boxes=%d\n",
	  state->verbose, state->alpha, state->mode, bins, boxes);
  fprintf(state->ostream, "\n       single.......iteration                   ");
  fprintf(state->ostream, "accumulated......results   \n");

  fprintf(state->ostream, "iteration     integral    sigma             integral   ");
  fprintf(state->ostream, "      sigma     chi-sq/it\n number\n");
  fflush(state->ostream);

}

void print_res(gsl_monte_vegas_state* state, 
	       int itr, double res, double err, double cum_res, double cum_err, 
	       double chi_sq)
{
  if (state->ostream == (FILE*) NULL) 
    return;

  fprintf(state->ostream, 
	  "%4d        %6.4e %10.2e          %6.4e      %8.2e  %10.2e\n",
	  itr, res, err, cum_res, cum_err, chi_sq);
  fflush(state->ostream);
}

void print_grid(gsl_monte_vegas_state* state, unsigned long dim)
{
  int mod, i, j;
  int p = state->verbose;
  if (p < 1 ) 
    return;
  if (state->ostream == (FILE*) NULL) 
    return;

  for (j = 0; j < dim; ++j) {
    fprintf(state->ostream, "\n axis %d \n", j);
    fprintf(state->ostream, "      x          delta         x     ");
    fprintf(state->ostream, "    delta         x          delta\n");
    for (i = 1 + p / 2, mod = 1; i <= state->bins; i += p, ++mod) {
      fprintf(state->ostream, "%11.2e%13.2e ", 
	      state->y_bin[i][j], state->bin_sum[i][j]);
      if (mod % 3 == 0)
	fprintf(state->ostream, "\n");
    }
    fprintf(state->ostream, "\n");
  }
  fprintf(state->ostream, "\n");
  fflush(state->ostream);

}
