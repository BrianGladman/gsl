#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <gsl_ran.h>
#include <gsl_siman.h>

/* A cute example of this program to then look at stuff would be:
  ./siman_test D3 | grep -v \# | awk '{print $3}' | sed 's/::/    /gp' > output

  then go into gnuplot and try typing:
  gnuplot> set parametric
  gnuplot> splot 'output' with lines
  */


/* my first cut at a simulated annealing optimizer */
void siman_solve(Element *x0_p, double (*Efunc)(Element x),
		 void (*take_step)(Element *x_p, double step_size),
		 double distance(Element x, Element y),
		 void print_position(Element x),
		 Ssiman_params params)
{
  /* the new set of trial points, and their energies and  probabilities */
  Element x, *new_x;
  double *energies, *probs, *sum_probs;
  double Ex;			/* energy of the chosen point */
  double T;			/* the temperature */
  int i, done;
  double throw;			/* throw the die to choose a new "x" */
  int n_iter;

  printf("#  iter    temperature       position         delta_pos        energy\n");

  new_x = (Element *) malloc(params.n_tries*sizeof(Element));
  energies = (double *) malloc(params.n_tries*sizeof(double));
  probs = (double *) malloc(params.n_tries*sizeof(double));
  sum_probs = (double *) malloc(params.n_tries*sizeof(double));

  gsl_ran_seed(time(0L));

  T = params.t_initial;
  x = *x0_p;
  done = 0;

  n_iter = 0;
  while (!done) {
    Ex = Efunc(x);
    for (i = 0; i < params.n_tries-1; ++i) { /* only go to N_TRIES-2 */
      /* center the new_x[] around x, then pass it to take_step() */
      sum_probs[i] = 0;
      new_x[i] = x;
      take_step(&(new_x[i]), params.step_size);
      energies[i] = (*Efunc)(new_x[i]);
      probs[i] = exp(-(energies[i]-Ex)/(params.k*T));
    }
    /* now add in the old value of "x", so it is a contendor */
    new_x[params.n_tries-1] = x;
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
	x = new_x[i];
	break;
      }
    }
    printf("%5d\t%12g\t", n_iter, T);
    print_position(x);
    printf("\t%12g\t%12g\n", distance(x, *x0_p), Ex);
    T /= params.mu_t;
    ++n_iter;
    if (T < params.t_min) {
      done = 1;
    }
  }
  /* now return the value via x0_p */
  *x0_p = x;
/*  printf("the result is: %g (E=%g)\n", x, Ex); */
}
