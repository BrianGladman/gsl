static void
compute_moments (double par, double * cheb);

static int
dgtsl (size_t n, double *c, double *d, double *e, double *b);

struct fn_fourier_params
{
  gsl_function *function;
  double omega;
};

static double fn_sin (double t, void *params);
static double fn_cos (double t, void *params);

static void compute_moments (double par, double *moment);

static void
qc25f (gsl_function * f, double a, double b, 
       gsl_integration_qawo_workspace * wf, size_t level,
       double *result, double *abserr, double *resabs, double *resasc);

static void
qc25f (gsl_function * f, double a, double b, 
       gsl_integration_qawo_workspace * wf, size_t level,
       double *result, double *abserr, double *resabs, double *resasc)
{
  const double center = 0.5 * (a + b);
  const double half_length = 0.5 * (b - a);
  const double omega = wf->omega ;
  
  const double par = omega * half_length;

  if (fabs (par) < 2)
    {
      gsl_function weighted_function;
      struct fn_fourier_params fn_params;

      fn_params.function = f;
      fn_params.omega = omega;

      if (wf->sine == GSL_INTEG_SINE) 
	{
	  weighted_function.function = &fn_sin;
	}
      else
	{
	  weighted_function.function = &fn_cos;
	}

      weighted_function.params = &fn_params;

      gsl_integration_qk15 (&weighted_function, a, b, result, abserr,
			    resabs, resasc);
      
      return;
    }
  else
    {
      double *moment;
      double cheb12[13], cheb24[25], chebmo[25];
      double result_abs, res12_cos, res12_sin, res24_cos, res24_sin;
      double est_cos, est_sin;
      double c, s;
      size_t i;

      gsl_integration_qcheb (f, a, b, cheb12, cheb24);

      if (level < wf->i)
	{
	  /* already in the cache */

	  moment = wf->chebmo + 25 * level;
	}
      else
	{
	  compute_moments (par, chebmo);
	  
	  moment = chebmo;

	  /* check for overflowing the table here (size n) */

	  if (level < wf->n)
	    {
	      /* fill in any intervening moments */

	      for (i = wf->i ; i < level; i++)
		{
		  compute_moments (par, wf->chebmo + 25*i);
		}
	      
	      /* copy the moments for the current level into the cache too */

	      for (i = 0; i < 25; i++)
		{
		  wf->chebmo[25 * level + i] = moment[i];
		}
	    }
	}


      res12_cos = cheb12[12] * moment[12];
      res12_sin = 0 ;

      for (i = 0; i < 6; i++)
	{
	  size_t k = 10 - 2 * i;
	  res12_cos += cheb12[k] * moment[k];
	  res12_sin += cheb12[k + 1] * moment[k + 1];
	}

      res24_cos = cheb24[24] * moment[24];
      res24_sin = 0 ;

      result_abs = fabs(cheb24[24]) ;

      for (i = 0; i < 12; i++)
	{
	  size_t k = 22 - 2 * i;
	  res24_cos += cheb24[k] * moment[k];
	  res24_sin += cheb24[k + 1] * moment[k + 1];
	  result_abs += fabs(cheb24[k]) + fabs(cheb24[k+1]);
	}

      est_cos = fabs(res24_cos - res12_cos);
      est_sin = fabs(res24_sin - res12_sin);

      c = half_length * cos(center * omega);
      s = half_length * sin(center * omega);

      if (wf->sine == GSL_INTEG_SINE)
	{
	  *result = c * res24_sin + s * res24_cos;
	  *abserr = fabs(c * est_sin) + fabs(s * est_cos);
	}
      else
	{
      	  *result = c * res24_cos - s * res24_sin;
	  *abserr = fabs(c * est_cos) + fabs(s * est_sin);
	}
      

      *resabs = result_abs * half_length;
      *resasc = GSL_DBL_MAX;

      return;
    }
}

static double
fn_sin (double x, void *params)
{
  struct fn_fourier_params *p = (struct fn_fourier_params *) params;
  gsl_function *f = p->function;
  double w = p->omega;
  double wx = w * x;
  double sinwx = sin(wx) ;
  return GSL_FN_EVAL (f, x) * sinwx;
}

static double
fn_cos (double x, void *params)
{
  struct fn_fourier_params *p = (struct fn_fourier_params *) params;
  gsl_function *f = p->function;
  double w = p->omega;
  double wx = w * x;
  double coswx = cos(wx) ;
  return GSL_FN_EVAL (f, x) * coswx ;
}

static void
compute_moments (double par, double *chebmo)
{
  double v[28], d[25], d1[25], d2[25];

  const size_t noeq = 25;
  
  const double par2 = par * par;
  const double par4 = par2 * par2;
  const double par22 = par2 + 2.0;

  const double sinpar = sin (par);
  const double cospar = cos (par);

  size_t i;

  /* compute the chebyschev moments with respect to cosine */

  double ac = 8 * cospar;
  double as = 24 * par * sinpar;

  v[0] = 2 * sinpar / par;
  v[1] = (8 * cospar + (2 * par2 - 8) * sinpar / par) / par2;
  v[2] = (32 * (par2 - 12) * cospar
	  + (2 * ((par2 - 80) * par2 + 192) * sinpar) / par) / par4;

  if (fabs (par) <= 24)
    {
      /* compute the moments as the solution of a boundary value
         problem using the asyptotic expansion as an endpoint */
      
      double an2, ass, asap;
      double an = 6;
      size_t k;

      for (k = 0; k < noeq - 1; k++)
	{
	  an2 = an * an;
	  d[k] = -2 * (an2 - 4) * (par22 - 2 * an2);
	  d2[k] = (an - 1) * (an - 2) * par2;
	  d1[k + 1] = (an + 3) * (an + 4) * par2;
	  v[k + 3] = as - (an2 - 4) * ac;
	  an = an + 2.0;
	}

      an2 = an * an;

      d[noeq - 1] = -2 * (an2 - 4) * (par22 - 2 * an2);
      v[noeq + 2] = as - (an2 - 4) * ac;
      v[3] = v[3] - 56 * par2 * v[2];

      ass = par * sinpar;
      asap = (((((210 * par2 - 1) * cospar - (105 * par2 - 63) * ass) / an2
		- (1 - 15 * par2) * cospar + 15 * ass) / an2 
	       - cospar + 3 * ass) / an2 
	      - cospar) / an2;
      v[noeq + 2] = v[noeq + 2] - 2 * asap * par2 * (an - 1) * (an - 2);

      dgtsl (noeq, d1, d, d2, v + 3);

    }
  else
    {
      /* compute the moments by forward recursion */
      size_t k;
      double an = 4;

      for (k = 3; k < 13; k++)
	{
	  double an2 = an * an;
	  v[k] = ((an2 - 4) * (2 * (par22 - 2 * an2) * v[k - 1] - ac)
		  + as - par2 * (an + 1) * (an + 2) * v[k - 2]) 
	    / (par2 * (an - 1) * (an - 2));
	  an = an + 2.0;
	}
    }


  for (i = 0; i < 13; i++)
    {
      chebmo[2 * i] = v[i];
    }

  /* compute the chebyschev moments with respect to sine */

  v[0] = 2 * (sinpar - par * cospar) / par2;
  v[1] = (18 - 48 / par2) * sinpar / par2 + (-2 + 48 / par2) * cospar / par;

  ac = -24 * par * cospar;
  as = -8 * sinpar;

  if (fabs (par) <= 24)
    {
      /* compute the moments as the solution of a boundary value
         problem using the asyptotic expansion as an endpoint */

      size_t k;
      double an2, ass, asap;
      double an = 5;

      for (k = 0; k < noeq - 1; k++)
	{
	  an2 = an * an;
	  d[k] = -2 * (an2 - 4) * (par22 - 2 * an2);
	  d2[k] = (an - 1) * (an - 2) * par2;
	  d1[k + 1] = (an + 3) * (an + 4) * par2;
	  v[k + 2] = ac + (an2 - 4) * as;
	  an = an + 2.0;
	}
      
      an2 = an * an;

      d[noeq - 1] = -2 * (an2 - 4) * (par22 - 2 * an2);
      v[noeq + 1] = ac + (an2 - 4) * as;
      v[2] = v[2] - 42 * par2 * v[1];

      ass = par * cospar;
      asap = (((((105 * par2 - 63) * ass - (210 * par2 - 1) * sinpar) / an2
		+ (15 * par2 - 1) * sinpar
		- 15 * ass) / an2 - sinpar - 3 * ass) / an2 - sinpar) / an2;
      v[noeq + 1] = v[noeq + 1] - 2 * asap * par2 * (an - 1) * (an - 2);

      dgtsl (noeq, d1, d, d2, v + 2);

    }
  else
    {
      /* compute the moments by forward recursion */
      size_t k;
      double an = 3;
      for (k = 2; k < 12; k++)
	{
	  double an2 = an * an;
	  v[k] = ((an2 - 4) * (2 * (par22 - 2 * an2) * v[k - 1] + as)
		  + ac - par2 * (an + 1) * (an + 2) * v[k - 2]) 
	    / (par2 * (an - 1) * (an - 2));
	  an = an + 2.0;
	}
    }

  for (i = 0; i < 12; i++)
    {
      chebmo[2 * i + 1] = v[i];
    }

}

static int
dgtsl (size_t n, double *c, double *d, double *e, double *b)
{
  /* solves a tridiagonal matrix A x = b 
     
     c[1 .. n - 1]   subdiagonal of the matrix A
     d[0 .. n - 1]   diagonal of the matrix A
     e[0 .. n - 2]   superdiagonal of the matrix A

     b[0 .. n - 1]   right hand side, replaced by the solution vector x */

  size_t k;

  c[0] = d[0];

  if (n == 0)
    {
      return GSL_SUCCESS;
    }

  if (n == 1)
    {
      b[0] = b[0] / d[0] ;
      return GSL_SUCCESS;
    }

  d[0] = e[0];
  e[0] = 0;
  e[n - 1] = 0;

  for (k = 0; k < n - 1; k++)
    {
      size_t k1 = k + 1;

      if (fabs (c[k1]) >= fabs (c[k]))
	{
	  {
	    double t = c[k1];
	    c[k1] = c[k];
	    c[k] = t;
	  };
	  {
	    double t = d[k1];
	    d[k1] = d[k];
	    d[k] = t;
	  };
	  {
	    double t = e[k1];
	    e[k1] = e[k];
	    e[k] = t;
	  };
	  {
	    double t = b[k1];
	    b[k1] = b[k];
	    b[k] = t;
	  };
	}

      if (c[k] == 0)
	{
	  return GSL_FAILURE ;
	}

      {
	double t = -c[k1] / c[k];

	c[k1] = d[k1] + t * d[k];
	d[k1] = e[k1] + t * e[k];
	e[k1] = 0;
	b[k1] = b[k1] + t * b[k];
      }

    }

  if (c[n - 1] == 0)
    {
      return GSL_FAILURE;
    }


  b[n - 1] = b[n - 1] / c[n - 1];

  b[n - 2] = (b[n - 2] - d[n - 2] * b[n - 1]) / c[n - 2];

  for (k = n ; k > 2; k--)
    {
      size_t kb = k - 3;
      b[kb] = (b[kb] - d[kb] * b[kb + 1] - e[kb] * b[kb + 2]) / c[kb];
    }

  return GSL_SUCCESS;
}
