#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gsl_rng.h>
#include <gsl_siman.h>

/* implementation of a basic simulated annealing algorithm */

void 
gsl_siman_solve (const gsl_rng * r, void *x0_p, gsl_Efunc_t Ef,
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
  int n_evals = 0, n_iter = 0, n_accepts, n_rejects, n_eless;

  distance = 0 ; /* This parameter is not currently used */

  x = (void *) malloc (element_size);
  new_x = (void *) malloc (element_size);

  T = params.t_initial;
  memcpy (x, x0_p, element_size);
  done = 0;

  if (print_position)
    {
      printf ("#-iter  #-evals   temperature     position   energy\n");
    }

  while (!done)
    {
      E = Ef (x);
      n_accepts = 0;
      n_rejects = 0;
      n_eless = 0;
      for (i = 0; i < params.n_tries - 1; ++i)
	{
	  memcpy (new_x, x, element_size);
	  take_step (r, new_x, params.step_size);
	  new_E = Ef (new_x);
	  ++n_evals;		/* keep track of Ef() evaluations */
	  /* now take the crucial step: see if the new point is accepted
	     or not, as determined by the boltzman probability */
	  if (new_E < E)
	    {
	      /* yay! take a step */
	      memcpy (x, new_x, element_size);
	      E = new_E;
	      ++n_eless;
	    }
	  else if (exp (-(E - new_E) / (params.k * T)) * gsl_rng_uniform (r) < 0.5)
	    {
	      /* yay! take a step */
	      memcpy (x, new_x, element_size);
	      E = new_E;
	      ++n_accepts;
	    }
	  else
	    {
	      ++n_rejects;
	    }
	}

      if (print_position)
	{			/* see if we need to print stuff as we go */
/*       printf("%5d %12g %5d %3d %3d %3d", n_iter, T, n_evals, */
/*           100*n_eless/n_steps, 100*n_accepts/n_steps, */
/*           100*n_rejects/n_steps); */
	  printf ("%5d   %7d  %12g", n_iter, n_evals, T);
	  print_position (x);
	  printf ("  %12g\n", E);
	}

      /* apply the cooling schedule to the temperature */
      /* FIXME: I should also introduce a cooling schedule for the n_tries */
      T /= params.mu_t;
      ++n_iter;
      if (T < params.t_min)
	{
	  done = 1;
	}
    }

  /* at the end, copy the result onto the initial point, so we pass it
     back to the caller */
  memcpy (x0_p, x, element_size);

  free (x);
  free (new_x);
}

/* implementation of a simulated annealing algorithm with many tries */

void 
gsl_siman_solve_many (const gsl_rng * r, void *x0_p, gsl_Efunc_t Ef,
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
  double u;			/* throw the die to choose a new "x" */
  int n_iter;

  if (print_position)
    {
      printf ("#-iter    temperature       position");
      printf ("         delta_pos        energy\n");
    }

  x = (void *) malloc (params.n_tries * element_size);
  new_x = (void *) malloc (params.n_tries * element_size);
  energies = (double *) malloc (params.n_tries * sizeof (double));
  probs = (double *) malloc (params.n_tries * sizeof (double));
  sum_probs = (double *) malloc (params.n_tries * sizeof (double));

  T = params.t_initial;
  memcpy (x, x0_p, element_size);
  done = 0;

  n_iter = 0;
  while (!done)
    {
      Ex = Ef (x);
      for (i = 0; i < params.n_tries - 1; ++i)
	{			/* only go to N_TRIES-2 */
	  /* center the new_x[] around x, then pass it to take_step() */
	  sum_probs[i] = 0;
	  memcpy ((char *)new_x + i * element_size, x, element_size);
	  take_step (r, (char *)new_x + i * element_size, params.step_size);
	  energies[i] = Ef ((char *)new_x + i * element_size);
	  probs[i] = exp (-(energies[i] - Ex) / (params.k * T));
	}
      /* now add in the old value of "x", so it is a contendor */
      memcpy ((char *)new_x + (params.n_tries - 1) * element_size, x, element_size);
      energies[params.n_tries - 1] = Ex;
      probs[params.n_tries - 1] = exp (-(energies[i] - Ex) / (params.k * T));

      /* now throw biased die to see which new_x[i] we choose */
      sum_probs[0] = probs[0];
      for (i = 1; i < params.n_tries; ++i)
	{
	  sum_probs[i] = sum_probs[i - 1] + probs[i];
	}
      u = gsl_rng_uniform (r) * sum_probs[params.n_tries - 1];
      for (i = 0; i < params.n_tries; ++i)
	{
	  if (u < sum_probs[i])
	    {
	      memcpy (x, (char *)new_x + i * element_size, element_size);
	      break;
	    }
	}
      if (print_position)
	{
	  printf ("%5d\t%12g\t", n_iter, T);
	  print_position (x);
	  printf ("\t%12g\t%12g\n", distance (x, x0_p), Ex);
	}
      T /= params.mu_t;
      ++n_iter;
      if (T < params.t_min)
	{
	  done = 1;
	}
    }

  /* now return the value via x0_p */
  memcpy (x0_p, x, element_size);

  /*  printf("the result is: %g (E=%g)\n", x, Ex); */

  free (x);
  free (new_x);
  free (energies);
  free (probs);
  free (sum_probs);
}
