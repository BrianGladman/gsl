#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

/* include these with quotes instead of angle brackets so that we get the
   right version during development */
#include "gsl_ran.h"
#include "gsl_siman.h"

/* A cute example of this program to then look at stuff would be:
./siman_test D3 | grep -v "^#" | awk '{print $3}' | sed 's/::/    /gp' > output

  then go into gnuplot and try typing:
  gnuplot> set parametric
  gnuplot> splot 'output' with lines

  Or you can plot from a pipe:
  gnuplot> plot '<../siman/siman_test | grep -v "^#"' using 1:2 with lines
  gnuplot> plot '<../siman/siman_test | grep -v "^#"' using 1:3 with lines
  and so forth.
  */

/* implementation of a basic simulated annealing algorithm */
void gsl_siman_solve(void *x0_p, gsl_Efunc_t Ef,
		     gsl_siman_step_t take_step,
		     gsl_siman_metric_t distance,
		     gsl_siman_print_t print_position,
		     size_t element_size,
		     gsl_siman_params_t params)
{
  void *x, *new_x;
  double E, new_E;
  int i, done;
  double T;

  x = (void *) malloc(element_size);
  new_x = (void *) malloc(element_size);

  gsl_ran_seed(time(0L));

  T = params.t_initial;
  memcpy(x, x0_p, element_size);
  done = 0;

  while (!done) {
    E = Ef(x);
    for (i = 0; i < params.n_tries-1; ++i) {
      take_step(new_x, params.step_size);
      new_E = Ef(new_x);
      /* now take the crucial step: see if the new point is accepted
	 or not, as determined by the boltzman probability */
      if ((new_E < E)
	  || (exp(-(E - new_E)/(params.k*T)) * gsl_ran_uniform() < 0.5)) {
	/* yay! take a step */
	x = new_x;
      }
    }

    if (print_position) {	/* see if we need to print stuff as we go */
      printf("%d %t");
      print_position(x);
      print_position(new_x);
      printf("  %g  %g\n", E, new_E);
    }
    if (T < params.t_min) {
      done = 1;
    }
  }

  free(x);
  free(new_x);
}

/* implementation of a simulated annealing algorithm with many tries */
void gsl_siman_solve_many(void *x0_p, gsl_Efunc_t Ef,
			  gsl_siman_step_t take_step,
			  gsl_siman_metric_t distance,
			  gsl_siman_print_t print_position,
			  size_t element_size,
			  gsl_siman_params_t params)
{
  /* the new set of trial points, and their energies and probabilities */
  void *x, *new_x;
  double *energies, *probs, *sum_probs;
  double Ex;			/* energy of the chosen point */
  double T;			/* the temperature */
  int i, done;
  double throw;			/* throw the die to choose a new "x" */
  int n_iter;

  if (print_position) {
    printf("#  iter    temperature       position");
    printf("         delta_pos        energy\n");
  }

  x = (void *) malloc(params.n_tries*element_size);
  new_x = (void *) malloc(params.n_tries*element_size);
  energies = (double *) malloc(params.n_tries*sizeof(double));
  probs = (double *) malloc(params.n_tries*sizeof(double));
  sum_probs = (double *) malloc(params.n_tries*sizeof(double));

  gsl_ran_seed(time(0L));

  T = params.t_initial;
  memcpy(x, x0_p, element_size);
  done = 0;

  n_iter = 0;
  while (!done) {
    Ex = Ef(x);
    for (i = 0; i < params.n_tries-1; ++i) { /* only go to N_TRIES-2 */
      /* center the new_x[] around x, then pass it to take_step() */
      sum_probs[i] = 0;
      memcpy(new_x + i*element_size, x, element_size);
      take_step(new_x + i*element_size, params.step_size);
      energies[i] = Ef(new_x + i*element_size);
      probs[i] = exp(-(energies[i]-Ex)/(params.k*T));
    }
    /* now add in the old value of "x", so it is a contendor */
    memcpy(new_x + (params.n_tries-1)*element_size, x, element_size);
    energies[params.n_tries-1] = Ex;
    probs[params.n_tries-1] = exp(-(energies[i]-Ex)/(params.k*T));

    /* now throw biased die to see which new_x[i] we choose */
    sum_probs[0] = probs[0];
    for (i = 1; i < params.n_tries; ++i) {
      sum_probs[i] = sum_probs[i-1] + probs[i];
    }
    throw = gsl_ran_uniform() * sum_probs[params.n_tries-1];
    for (i = 0; i < params.n_tries; ++i) {
      if (throw < sum_probs[i]) {
	memcpy(x, new_x + i*element_size, element_size);
	break;
      }
    }
    if (print_position) {
      printf("%5d\t%12g\t", n_iter, T);
      print_position(x);
      printf("\t%12g\t%12g\n", distance(x, x0_p), Ex);
    }
    T /= params.mu_t;
    ++n_iter;
    if (T < params.t_min) {
      done = 1;
    }
  }
  /* now return the value via x0_p */
  memcpy(x0_p, x, element_size);
/*  printf("the result is: %g (E=%g)\n", x, Ex); */
  free(x);
  free(new_x);
  free(energies);
  free(probs);
  free(sum_probs);
}
