static int
set (void *vstate, gsl_multifit_function_fdf * fdf, gsl_vector * x,
     gsl_vector * f, gsl_vector * dx, const gsl_vector * weights,
     int scale)
{
  lmder_state_t *state = (lmder_state_t *) vstate;

  gsl_matrix *r = state->r;
  gsl_vector *tau = state->tau;
  gsl_vector *qtf = state->qtf;
  gsl_vector *diag = state->diag;
  gsl_vector *work1 = state->work1;
  gsl_permutation *perm = state->perm;

  int signum;

  /* start counting function and Jacobian evaluations */
  fdf->nevalf = 0;
  fdf->nevaldf = 0;

  /* keep pointer to weight vector */
  state->weights = weights;

  /* return immediately if evaluation raised error */
  {
    int status;

    /* Evaluate function at x */
    status = GSL_MULTIFIT_FN_EVAL_F (fdf, x, f);
    if (status)
      return status;

    /* Evaluate Jacobian at x and store in state->r */
    if (fdf->df)
      status = GSL_MULTIFIT_FN_EVAL_DF (fdf, x, r);
    else /* finite difference approximation */
      status = gsl_multifit_fdfsolver_dif_df(x, fdf, f, r);

    if (status)
      return status;
  }

  state->par = 0;
  state->iter = 1;
  state->fnorm = enorm (f);

  gsl_vector_set_all (dx, 0.0);

  /* store column norms in diag */

  if (scale)
    {
      compute_diag (r, diag);
    }
  else
    {
      gsl_vector_set_all (diag, 1.0);
    }

  /* set delta to 100 |D x| or to 100 if |D x| is zero */

  state->xnorm = scaled_enorm (diag, x);
  state->delta = compute_delta (diag, x);

  /* Factorize J = Q R P^T */
  gsl_linalg_QRPT_decomp (r, tau, perm, &signum, work1);

  /* compute qtf = Q^T f */
  gsl_vector_memcpy (qtf, f);
  gsl_linalg_QR_QTvec (r, tau, qtf);

  gsl_vector_set_zero (state->rptdx);
  gsl_vector_set_zero (state->w);

  /* Zero the trial vector, as in the alloc function */

  gsl_vector_set_zero (state->f_trial);

#ifdef DEBUG
  printf("r = "); gsl_matrix_fprintf(stdout, r, "%g");
  printf("perm = "); gsl_permutation_fprintf(stdout, perm, "%d");
  printf("tau = "); gsl_vector_fprintf(stdout, tau, "%g");
#endif

  return GSL_SUCCESS;
}
