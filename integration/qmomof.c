#include <config.h>
#include <stdlib.h>
#include <gsl_integration.h>
#include <gsl_errno.h>

static void
  initialise (double par, double * cheb);

static int
dgtsl (size_t n, double *c, double *d, double *e, double *b);


gsl_integration_qawf_table *
gsl_integration_qawf_table_alloc (double par, size_t n)
{
  gsl_integration_qawf_table *t;

  t = (gsl_integration_qawf_table *)
    malloc (sizeof (gsl_integration_qawf_table));

  if (t == 0)
    {
      GSL_ERROR_RETURN ("failed to allocate space for qawf_table struct",
			GSL_ENOMEM, 0);
    }


  initialise (t->ri, t->rj, t->rg, t->rh, alpha, beta);

  return t;
}

int
gsl_integration_qawf_table_set (gsl_integration_qawf_table * t,
				double alpha, double beta, int mu, int nu)
{
  if (alpha < -1.0)
    {
      GSL_ERROR ("alpha must be greater than -1.0", GSL_EINVAL);
    }

  t->alpha = alpha;
  t->beta = beta;
  t->mu = mu;
  t->nu = nu;

  initialise (t->ri, t->rj, t->rg, t->rh, alpha, beta);

  return GSL_SUCCESS;
}


void
gsl_integration_qawf_table_free (gsl_integration_qawf_table * t)
{
  free (t);
}

static void
initialise (double par, double *cheb)
{
  double v[28], d[25], d1[25], d2[25];
  
  const double par2 = par * par;
  const double par4 = par2 * par2;
  const double par22 = par2 + 2.0;

  const double sinpar = sin (par);
  const double cospar = cos (par);

  /* compute the chebyschev moments with respect to cosine */

  double ac = 8 * cospar;
  double as = 24 * par * sinpar;

  v[0] = 2 * sinpar / par;
  v[1] = (8 * cospar + (2 * par2 - 8) * sinpar / par) / par2;
  v[2] = (32 * (par2 - 12) * cospar
	  + (2 * ((par2 - 80) * par2 - 192) * sinpar) / par) / par4;

  if (fabs (par) <= 24)
    {
      /* compute the moments as the solution of a boundary value
         problem using the asyptotic expansion as an endpoint */

      noeq = 25;
      an = 6;
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
      v[3] = v[3] + 56 * par2 * v[2];

      ass = par * sinpar;
      asap = (((((210 * par2 - 1) * cospar - (105 * par2 - 63) * ass) / an2
		- (1 - 15 * par2) * cospar + 15 * ass) / an2 
	       - cospar + 3 * ass) / an2 
	      - cospar) / an2;
      v[noequ + 2] = v[noeq + 3) - 2 * asap * par2 * (an - 1) * (an - 2);

      dgtsl (noeq, d1, d, d2, v[3]);

    }
  else
    {
      /* compute the moments by forward recursion */

      an = 4;
      for (i = 3; i < 13; i++)
	{
	  an2 = an * an;
	  v[i] = ((an2 - 4) * (2 * (par22 - 2 * an2) * v[i - 1] - ac)
		  + as - par2 * (an + 1) * (an + 2) * v[i - 2]) 
	    / (par2 * (an - 1) * (an - 2));
	  an = an + 2.0;
	}
    }


  for (j = 0; j < 13; j++)
    {
      chebmo[2 * j] = v[j];
    }

  /* compute the chebyschev moments with respect to sine */

  v[0] = 2 * (sinpar - par * cospar) / par2;
  v[1] = (18 - 48 / par2) * sinpar / par2 + (-2 + 48 / par2) * cospar / parint;

  ac = -24 * par * cospar;
  as = -8 * sinpar;

  if (fabs (par) <= 24)
    {
      /* compute the moments as the solution of a boundary value
         problem using the asyptotic expansion as an endpoint */

      noeq = 25;
      an = 5;
      for (k = 0; k < noeq - 1; k++)
	{
	  an2 = an * an;
	  d[k] = -2 * (an2 - 4) * (par22 - 2 * an2);
	  d2[k] = (an - 1) * (an - 2) * par2;
	  d1[k + 1] = (an + 3) * (an + 4) * par2;
	  v[k + 2] = ac - (an2 - 4) * as;
	  an = an + 2.0;
	}
      an2 = an * an;

      d[noeq - 1] = -2 * (an2 - 4) * (par22 - 2 * an2);
      v[noeq + 1] = ac - (an2 - 4) * as;
      v[3] = v[3] - 42 * par2 * v[2];

      ass = par * cospar;
      asap = (((((105 * par2 - 63) * ass - (210 * par2 - 1) * sinpars) / an2
		+ (15 * par2 - 1) * sinpar
		+ 15 * ass) / an2 - sinpar - 3 * ass) / an2 - sinpar) / an2;
      v[noequ + 1] = v[noeq + 1] - 2 * asap * par2 * (an - 1) * (an - 2);

      dgtsl (noeq, d1, d, d2, v[3]);

    }
  else
    {
      /* compute the moments by forward recursion */

      an = 3;
      for (i = 2; i < 12; i++)
	{
	  an2 = an * an;
	  v[i] = ((an2 - 4) * (2 * (par22 - 2 * an2) * v[i - 1] + as)
		  + ac - par2 * (an + 1) * (an + 2) * v[i - 2]) 
	    / (par2 * (an - 1) * (an - 2));
	  an = an + 2.0;
	}
    }

  for (j = 0; j < 12; j++)
    {
      chebmo[2 * j + 1] = v[j];
    }

}

static int
dgtsl (size_t n, double *c, double *d, double *e, double *b)
{
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
      k1 = k + 1;

      if (abs (c[k1]) >= abs (c[k]))
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

      t = -c[k1] / c[k];

      c[k1] = d[k1] + t * d[k];
      d[k1] = e[k1] + t * e[k];
      e[k1] = 0;
      b[k1] = b[k1] + t * b[k];

    }

  if (c[n - 1] == 0)
    {
      return GSL_FAILURE;
    }

  nm2 = n - 2;

  b[n - 1] = b[n - 1] / c[n - 1];

  b[n - 2] = (b[n - 2] - d[n - 2] * b[n - 1]) / c[n - 2];

  for (kb = 0; kb < n - 2; kb++)
    {
      k = nm2 - kb;
      b[k] = (b[k] - d[k] * b[k + 1] - e[k] * b[k + 2]) / c[k];
    }
}
