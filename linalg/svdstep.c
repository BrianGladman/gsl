inline static void
chop_small_elements (gsl_vector * d, gsl_vector * f)
{
  const size_t N = d->size;
  double d_i = gsl_vector_get (d, 0);

  size_t i;

  for (i = 0; i < N - 1; i++)
    {
      double f_i = gsl_vector_get (f, i);
      double d_ip1 = gsl_vector_get (d, i + 1);

      if (fabs (f_i) < GSL_DBL_EPSILON * (fabs (d_i) + fabs (d_ip1)))
	{
	  gsl_vector_set (f, i, 0.0);
	}
      d_i = d_ip1;
    }

}

#ifdef JUNK
inline static void
chase_out_zero (gsl_vector * d, gsl_vector * f, double gc[], double gs[])
{
  const size_t n = d->size;

  double c = 0.0, s = 1.0;

  size_t i;

  for (i = 0; i < N - 1; i++)
    {
      double d_i = gsl_vector_get (d, i);
      double f_i = gsl_vector_get (f, i);

      double f1 = s * f_i;

      gsl_vector_set (f, i, c * f_i);

      create_givens (d_i, f1, &c, &s);

      gsl_vector_set (d, i, hypot (d_i, f1));

      if (gc)
	gc[i] = c;
      if (gs)
	gs[i] = s;
    }
}
#endif

static double
trailing_eigenvalue (const gsl_vector * d, const gsl_vector * f)
{
  const size_t n = d->size;

  double da = gsl_vector_get (d, n - 2);
  double db = gsl_vector_get (d, n - 1);
  double fa = (n > 2) ? gsl_vector_get (f, n - 3) : 0.0;
  double fb = gsl_vector_get (f, n - 2);

  double ta = da * da + fa * fa;
  double tb = db * db + fb * fb;
  double tab = da * fb;

  double dt = (ta - tb) / 2.0;

  double mu;

  if (dt > 0)
    {
      mu = tb - (tab * tab) / (dt + hypot (dt, tab));
    }
  else if (dt < 0)
    {
      mu = tb + (tab * tab) / ((-dt) + hypot (dt, tab));
    }
  else
    {
      /* FIXME: what is the best choice of mu for dt == 0 ? */
      mu = tb;
    }

  return mu;
}

static void
qrstep (gsl_vector * d, gsl_vector * f, gsl_matrix * U, gsl_matrix * V)
{
  const size_t M = U->size1;
  const size_t N = V->size1;
  const size_t n = d->size;
  double y, z;
  double ak, bk, zk, ap, bp, aq, bq;
  size_t i, k;

  double d0 = gsl_vector_get (d, 0);
  double f0 = gsl_vector_get (f, 0);

  double d1 = gsl_vector_get (d, 1);
  double f1 = (n > 2) ? gsl_vector_get (f, 1) : 0.0;

  double mu = trailing_eigenvalue (d, f);

  y = d0 * d0 - mu;
  z = d0 * f0;

  /* Set up the recurrence for Givens rotations on a bidiagonal matrix */

  ak = 0;
  bk = 0;

  ap = d0;
  bp = f0;

  aq = d1;
  bq = f1;

  for (k = 0; k < n - 1; k++)
    {
      double c, s;
      create_givens (y, z, &c, &s);

      /* Compute V <= V G */

      for (i = 0; i < N; i++)
	{
	  double Vip = gsl_matrix_get (V, i, k);
	  double Viq = gsl_matrix_get (V, i, k + 1);
	  gsl_matrix_set (V, i, k, c * Vip - s * Viq);
	  gsl_matrix_set (V, i, k + 1, s * Vip + c * Viq);
	}

      /* compute B <= B G */

      {
	double bk1 = c * bk - s * z;

	double ap1 = c * ap - s * bp;
	double bp1 = s * ap + c * bp;
	double zp1 = -s * aq;

	double aq1 = c * aq;

	if (k > 0)
	  {
	    gsl_vector_set (f, k - 1, bk1);
	  }

	ak = ap1;
	bk = bp1;
	zk = zp1;

	ap = aq1;

	if (k < n - 2)
	  {
	    bp = gsl_vector_get (f, k + 1);
	  }
	else
	  {
	    bp = 0.0;
	  }

	y = ak;
	z = zk;
      }

      create_givens (y, z, &c, &s);

      /* Compute U <= U G */

      for (i = 0; i < M; i++)
	{
	  double Uip = gsl_matrix_get (U, i, k);
	  double Uiq = gsl_matrix_get (U, i, k + 1);
	  gsl_matrix_set (U, i, k, c * Uip - s * Uiq);
	  gsl_matrix_set (U, i, k + 1, s * Uip + c * Uiq);
	}

      /* compute B <= G^T B */

      {
	double ak1 = c * ak - s * zk;
	double bk1 = c * bk - s * ap;
	double zk1 = -s * bp;

	double ap1 = s * bk + c * ap;
	double bp1 = c * bp;

	gsl_vector_set (d, k, ak1);

	ak = ak1;
	bk = bk1;
	zk = zk1;

	ap = ap1;
	bp = bp1;

	if (k < n - 2)
	  {
	    aq = gsl_vector_get (d, k + 2);
	  }
	else
	  {
	    aq = 0.0;
	  }

	y = bk;
	z = zk;
      }
    }

  gsl_vector_set (f, n - 2, bk);
  gsl_vector_set (d, n - 1, ap);
}
