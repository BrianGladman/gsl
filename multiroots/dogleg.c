#include "enorm.c"

static int
newton_direction (const gsl_matrix * r, const gsl_vector * qtf, gsl_vector * p)
{
  const size_t N = r->size2;
  size_t i;
  int status;

  gsl_vector_copy (p, qtf);

  status = gsl_la_Rsolve_impl (r, p);

  for (i = 0; i < N; i++)
    {
      double pi = gsl_vector_get (p, i);
      gsl_vector_set (p, i, -pi);
    }

  return status;
}

static void
gradient_direction (const gsl_matrix * r, const gsl_vector * qtf,
		    const gsl_vector * diag, gsl_vector * g)
{
  const size_t M = r->size1;
  const size_t N = r->size2;

  size_t i, j;

  for (j = 0; j < M; j++)
    {
      double sum = 0;
      double dj;

      for (i = 0; i < N; i++)
	{
	  sum += gsl_matrix_get (r, i, j) * gsl_vector_get (qtf, i);
	}

      dj = gsl_vector_get (diag, j);
      gsl_vector_set (g, j, -sum / dj);
    }
}

static void
minimum_step (double gnorm, const gsl_vector * diag, gsl_vector * g)
{
  const size_t N = g->size;
  size_t i;

  for (i = 0; i < N; i++)
    {
      double gi = gsl_vector_get (g, i);
      double di = gsl_vector_get (diag, i);
      gsl_vector_set (g, i, gi / (di * gnorm));
    }
}

static void
compute_Rg (const gsl_matrix * r, const gsl_vector * gradient, gsl_vector * Rg)
{
  const size_t N = r->size2;
  size_t i, j;

  for (i = 0; i < N; i++)
    {
      double sum = 0;

      for (j = i; j < N; j++)
	{
	  double gj = gsl_vector_get (gradient, j);
	  double rij = gsl_matrix_get (r, i, j);
	  sum += rij * gj;
	}

      gsl_vector_set (Rg, i, sum);
    }
}

static void
scaled_addition (double alpha, gsl_vector * newton, double beta, gsl_vector * gradient, gsl_vector * p)
{
  const size_t N = p->size;
  size_t i;

  for (i = 0; i < N; i++)
    {
      double ni = gsl_vector_get (newton, i);
      double gi = gsl_vector_get (gradient, i);
      gsl_vector_set (p, i, alpha * ni + beta * gi);
    }
}

static int
dogleg (const gsl_matrix * r, const gsl_vector * qtf,
	const gsl_vector * diag, double delta,
	gsl_vector * newton, gsl_vector * gradient, gsl_vector * p)
{
  double qnorm, gnorm, sgnorm, bnorm, temp;

  newton_direction (r, qtf, newton);

  qnorm = scaled_enorm (diag, p);

  if (qnorm <= delta)
    {
      gsl_vector_copy (p, newton);
      return GSL_SUCCESS;
    }

  gradient_direction (r, qtf, diag, gradient);

  gnorm = enorm (gradient);

  if (gnorm == 0)
    {
      double alpha = delta / qnorm;
      double beta = 0;
      scaled_addition (alpha, newton, 0, gradient, p);
      return GSL_SUCCESS;
    }

  minimum_step (gnorm, diag, gradient);

  compute_Rg (r, gradient, p);	/* use p as temporary space to compute Rg */

  temp = enorm (p);
  sgnorm = (gnorm / temp) / temp;

  if (sgnorm > delta)
    {
      double alpha = 0;
      double beta = delta;
      scaled_addition (0, newton, beta, gradient, p);
      return GSL_SUCCESS;
    }

  bnorm = enorm (qtf);

  {
    double bg = bnorm / gnorm;
    double bq = bnorm / qnorm;
    double dq = delta / qnorm;
    double dq2 = dq * dq;
    double sd = sgnorm / delta;
    double sd2 = sd * sd;

    double t1 = bg * bq * sd;
    double t2 = t1 - dq * sd2 + sqrt ((t1 - dq2) * (1 - sd2));

    double alpha = dq * (1 - sd2) / t2;
    double beta = (1 - alpha) * sgnorm;

    scaled_addition (alpha, newton, beta, gradient, p);
  }

  return GSL_SUCCESS;
}
