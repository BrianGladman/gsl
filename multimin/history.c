
gsl_multimin_fdf_history *
gsl_multimin_fdf_history_alloc (size_t n)
{
  gsl_multimin_fdf_history *h;

  h = (gsl_multimin_fdf_history *) malloc (sizeof (gsl_multimin_fdf_history));

  if (h == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for multimin history struct",
		     GSL_ENOMEM, 0);
    }

  h->x = gsl_vector_calloc (n);

  if (h->x == 0)
    {
      free (h);
      GSL_ERROR_VAL ("failed to allocate space for x", GSL_ENOMEM, 0);
    }

  h->x1 = gsl_vector_calloc (n);

  if (h->x1 == 0)
    {
      free (h);
      gsl_vector_free (h->x);
      GSL_ERROR_VAL ("failed to allocate space for x1", GSL_ENOMEM, 0);
    }

  h->g = gsl_vector_calloc (n);

  if (h->g == 0)
    {
      free (h);
      gsl_vector_free (h->x);
      gsl_vector_free (h->x1);
      GSL_ERROR_VAL ("failed to allocate space for g", GSL_ENOMEM, 0);
    }

  h->g1 = gsl_vector_calloc (n);

  if (h->g1 == 0)
    {
      free (h);
      gsl_vector_free (h->x);
      gsl_vector_free (h->x1);
      gsl_vector_free (h->g);
      GSL_ERROR_VAL ("failed to allocate space for g1", GSL_ENOMEM, 0);
    }

  return h;
}

int
gsl_multimin_fdf_history_set (gsl_multimin_fdf_history * h,
			      gsl_multimin_function_fdf * fdf,
			      const gsl_vector * x)
{
  gsl_vector_memcpy (h->x, x);
  GSL_MULTIMIN_FN_EVAL_F_DF (fdf, h->x, &(h->f), h->g);
  if (!finite (h->f))
    GSL_ERROR ("function not continuous", GSL_EBADFUNC);
  return GSL_SUCCESS;
}

int
gsl_multimin_fdf_history_set_with_value (gsl_multimin_fdf_history * h,
					 gsl_multimin_function_fdf * fdf,
					 const gsl_vector * x, double fx)
{
  gsl_vector_memcpy (h->x, x);
  GSL_MULTIMIN_FN_EVAL_DF (fdf, h->x, h->g);
  h->f = fx;
  return GSL_SUCCESS;
}

int
gsl_multimin_fdf_history_step (gsl_multimin_fdf_history * h,
			       gsl_multimin_function_fdf * fdf,
			       const gsl_vector * direction, double step)
{
  gsl_vector_memcpy (h->g1, h->g);
  gsl_vector_memcpy (h->x1, h->x);
  h->f1 = h->f;
  gsl_multimin_compute_evaluation_point (h->x, h->x1, step, direction);
  GSL_MULTIMIN_FN_EVAL_F_DF (fdf, h->x, &(h->f), h->g);
  if (!finite (h->f))
    GSL_ERROR ("function not continuous", GSL_EBADFUNC);
  return GSL_SUCCESS;
}

int
gsl_multimin_fdf_history_step_with_value (gsl_multimin_fdf_history * h,
					  gsl_multimin_function_fdf * fdf,
					  const gsl_vector * direction,
					  double step, double fx)
{
  gsl_vector_memcpy (h->g1, h->g);
  gsl_vector_memcpy (h->x1, h->x);
  h->f1 = h->f;
  gsl_multimin_compute_evaluation_point (h->x, h->x1, step, direction);
  GSL_MULTIMIN_FN_EVAL_DF (fdf, h->x, h->g);
  h->f = fx;
  return GSL_SUCCESS;
}

void
gsl_multimin_fdf_history_free (gsl_multimin_fdf_history * h)
{
  gsl_vector_free (h->x);
  gsl_vector_free (h->x1);
  gsl_vector_free (h->g);
  gsl_vector_free (h->g1);
  free (h);
}
