static int
iterate (void *vstate, gsl_multifit_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx, int scale)
{
  lmder_state_t *state = (lmder_state_t *) vstate;

  const double xnorm = state->xnorm;
  const double fnorm = state->fnorm;

  gsl_matrix *q = state->q;
  gsl_matrix *r = state->r;
  gsl_vector *tau = state->tau;
  gsl_vector *diag = state->diag;
  gsl_vector *qtf = state->qtf;
  gsl_vector *x_trial = state->x_trial;
  gsl_vector *f_trial = state->f_trial;
  gsl_vector *df = state->df;
  gsl_vector *qtdf = state->qtdf;
  gsl_vector *rptdx = state->rptdx;
  gsl_vector *gradient = state->gradient;
  gsl_permutation *perm = state->perm;

  double prered, actred;
  double pnorm, fnorm1, fnorm1p, gnorm;
  double ratio;
  double dirder, delta, par;

  double p1 = 0.1, p25 = 0.25, p5 = 0.5, p75 = 0.75, p0001 = 0.0001;

  /* Compute qtf = Q^T f */

  compute_qtf (q, f, qtf);

  /* Compute norm of scaled gradient */

  gradient_direction (r, qtf, diag, perm, gradient);

  { 
    size_t iamax = gsl_blas_idamax (gradient);

    gnorm = fabs(gsl_vector_get (gradient, iamax) / fnorm);
  }

  /* Determine the levenberg-marquardt parameter */

lm_iteration:
  
  lmpar (r, qtf, diag, state->delta, state->newton, gradient, dx);

  /* Take a trial step */

  compute_trial_step (x, dx, state->x_trial);

  pnorm = scaled_enorm (diag, dx);

  if (state->iter == 1)
    {
      if (pnorm < state->delta)
	{
	  state->delta = pnorm;
	}
    }

  /* Evaluate function at x + p */

  GSL_MULTIFIT_FN_EVAL_F (fdf, x_trial, f_trial);

  fnorm1 = enorm (f_trial);

  /* Compute the scaled actual reduction */

  actred = compute_actual_reduction (fnorm, fnorm1);

  /* Compute rptdx = R P^T dx, noting that |J dx| = |R P^T dx| */

  compute_rptdx (r, perm, dx, rptdx);

  fnorm1p = enorm (rptdx);

  /* Compute the scaled predicted reduction = |J dx|^2 + 2 par |D dx|^2 */

  { 
    double t1 = fnorm1p / fnorm;
    double t2 = (sqrt(par) * pnorm) / fnorm;
    
    prered = t1 * t1 + t2 * t2 / p5;
    dirder = -(t1 * t1 + t2 * t2);
  }

  /* compute the ratio of the actual to predicted reduction */

  if (prered > 0)
    {
      ratio = actred / prered;
    }
  else
    {
      ratio = 0;
    }

  /* update the step bound */

  if (ratio > p25)
    {
      if (par == 0 || ratio >= p75)
        {
          delta = pnorm / p5;
          par *= p5;
        }
    }
  else
    {
      double temp = (actred < 0) ? p5*dirder / (dirder + p5 * actred) : p5;

      if (p1 * fnorm1 >= fnorm || temp < p1 ) temp = p1;

      par /= temp;
    }

  /* test for successful iteration, termination and stringent tolerances */

  if (ratio >= p0001)
    {
      gsl_vector_memcpy (x, x_trial);
      gsl_vector_memcpy (f, f_trial);

      GSL_MULTIFIT_FN_EVAL_DF (fdf, x_trial, J);

      /* wa2_j  = diag_j * x_j */
      state->xnorm = enorm(x);
      state->fnorm = fnorm1;
      state->iter++;

      /* Rescale if necessary */

      if (scale)
        {
          update_diag (J, diag);
        }
    }
  else if (fabs(actred) <= GSL_DBL_EPSILON  && prered <= GSL_DBL_EPSILON 
           && p5 * ratio <= 1.0)
    {
      return GSL_ETOLF ;
    }
  else if (delta <= GSL_DBL_EPSILON * state->xnorm)
    {
      return GSL_ETOLX;
    }
  else if (gnorm <= GSL_DBL_EPSILON)
    {
      return GSL_ETOLG;
    }
  else 
    {
      /* Repeat inner loop if unsuccessful */
      
      goto lm_iteration;
    }

  return GSL_SUCCESS;
}
