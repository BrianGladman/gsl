/*   Set of programs for printing vegas results   */
/* RCS $Id$ */

#include <stdio.h>
#include "gsl_vegas.h"
#include "gsl_vegas_print.h"
/* predeclare functions */

/* The following variables are needed for the print functions. */

FILE               *o_file;

void vegas_open_log(void)
{
  o_file = fopen("vegas.out", "a+");
}

void vegas_close_log(void)
{
  fclose(o_file);
}

void prn_lim(double a[], double b[], int m)
{
  int j;

  fprintf(o_file, "The limits of integration are:\n");
  for (j = 0; j < m; ++j)
    fprintf(o_file, "\nxl[%d]=%f    xu[%d]=%f", j, a[j], j, b[j]);
  fprintf(o_file, "\n");
  fflush(o_file);
}
void prn_head(gsl_monte_vegas_state* state, 
	      int num_dim, int calls, int it_num, int bins, int boxes)
{
  fprintf(o_file, "\nnum_dim=%d, calls=%d, it_num=%d, max_it_num=%d, acc=%.3f, ",
	  num_dim, calls, it_num, state->max_it_num, state->acc);
  fprintf(o_file, "verb=%d, alph=%.2f,\nmode=%d, bins=%d, boxes=%d\n",
	  state->verbose, state->alpha, state->mode, bins, boxes);
  fprintf(o_file, "\n       single.......iteration                   ");
  fprintf(o_file, "accumulated......results   \n");

  fprintf(o_file, "iteration     integral    sigma             integral   ");
  fprintf(o_file, "      sigma     chi-sq/it\n number\n");
  fflush(o_file);

}

void prn_res(int a, double b, double c, double d, double e, double f)
{
  fprintf(o_file, "%4d        %6.4e %10.2e          %6.4e      %8.2e  %10.2e\n",
	  a, b, c, d, e, f);
  fflush(o_file);
}

void prn_grid(double y[GSL_V_BINS_MAX/2][GSL_V_MAX_DIM], 
	      double s[GSL_V_BINS_MAX/2][GSL_V_MAX_DIM], 
	      int m, int n, int p)
{
  int mod, i, j;

  if (p < 1 ) 
    return;

  for (j = 0; j < m; ++j) {
    fprintf(o_file, "\n axis %d \n", j);
    fprintf(o_file, "      x          delta         x     ");
    fprintf(o_file, "    delta         x          delta\n");
    for (i = 1 + p / 2, mod = 1; i <= n; i += p, ++mod) {
      fprintf(o_file, "%11.2e%13.2e ", y[i][j], s[i][j]);
      if (mod % 3 == 0)
	fprintf(o_file, "\n");
    }
    fprintf(o_file, "\n");
  }
  fprintf(o_file, "\n");
  fflush(o_file);

}
