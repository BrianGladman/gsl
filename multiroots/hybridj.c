#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_multiroots.h>
#include <gsl_linalg.h>

#include "dogleg.c"

typedef struct
  {
    double delta;
    gsl_matrix * q;
    gsl_matrix * r;
    gsl_vector * rdiag;
    gsl_vector * diag;
    gsl_vector * qtf;
    gsl_vector * newton;
    gsl_vector * gradient;
  }
hybridj_state_t;

int hybridj_alloc (void * vstate, size_t n);
int hybridj_set (void * vstate, gsl_multiroot_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx);
int hybridj_iterate (void * vstate, gsl_multiroot_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx);
void hybridj_free (void * vstate);

int
hybridj_alloc (void * vstate, size_t n)
{
  hybridj_state_t * state = (hybridj_state_t *) vstate;
  gsl_matrix * q, *r;
  gsl_vector * rdiag, *diag, *qtf, *newton, *gradient;

  q = gsl_matrix_calloc (n,n);
  
  if (q == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for q", GSL_ENOMEM, 0);
    }

  state->q = q ;

  r = gsl_matrix_calloc (n,n);

  if (r == 0)
    {
      gsl_matrix_free(q);

      GSL_ERROR_RETURN ("failed to allocate space for r", GSL_ENOMEM, 0);
    }

  state->r = r ;

  rdiag = gsl_vector_calloc (n);

  if (rdiag == 0)
    {
      gsl_matrix_free(q);
      gsl_matrix_free(r);

      GSL_ERROR_RETURN ("failed to allocate space for rdiag", GSL_ENOMEM, 0);
    }

  state->rdiag = rdiag ;

  diag = gsl_vector_calloc (n);

  if (diag == 0)
    {
      gsl_matrix_free(q);
      gsl_matrix_free(r);
      gsl_vector_free(rdiag);

      GSL_ERROR_RETURN ("failed to allocate space for diag", GSL_ENOMEM, 0);
    }

  state->diag = diag ;

  qtf = gsl_vector_calloc (n);

  if (qtf == 0)
    {
      gsl_matrix_free(q);
      gsl_matrix_free(r);
      gsl_vector_free(rdiag);
      gsl_vector_free(diag);

      GSL_ERROR_RETURN ("failed to allocate space for qtf", GSL_ENOMEM, 0);
    }

  state->qtf = qtf ;

  newton = gsl_vector_calloc (n);

  if (newton == 0)
    {
      gsl_matrix_free(q);
      gsl_matrix_free(r);
      gsl_vector_free(rdiag);
      gsl_vector_free(diag);
      gsl_vector_free(qtf);

      GSL_ERROR_RETURN ("failed to allocate space for newton", GSL_ENOMEM, 0);
    }

  state->newton = newton ;

  gradient = gsl_vector_calloc (n);

  if (gradient == 0)
    {
      gsl_matrix_free(q);
      gsl_matrix_free(r);
      gsl_vector_free(rdiag);
      gsl_vector_free(diag);
      gsl_vector_free(qtf);
      gsl_vector_free(newton);

      GSL_ERROR_RETURN ("failed to allocate space for gradient", GSL_ENOMEM, 0);
    }

  state->gradient = gradient ;

  x_trial = gsl_vector_calloc (n);

  if (x_trial == 0)
    {
      gsl_matrix_free(q);
      gsl_matrix_free(r);
      gsl_vector_free(rdiag);
      gsl_vector_free(diag);
      gsl_vector_free(qtf);
      gsl_vector_free(newton);
      gsl_vector_free(gradient);

      GSL_ERROR_RETURN ("failed to allocate space for x_trial", GSL_ENOMEM, 0);
    }

  state->x_trial = x_trial;


  return GSL_SUCCESS;
}

int 
hybridj_set (void * vstate, gsl_multiroot_function_fdf * FDF, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx)
{
  hybridj_state_t * state = (hybridj_state_t *) vstate;

  size_t i, j, n = FDF->n ;

  gsl_matrix * q = state->q;
  gsl_matrix * r = state->r;
  gsl_vector * rdiag = state->rdiag;
  gsl_vector * diag = state->diag;
  gsl_vector * qtf = state->qtf;

  double Dx, factor;

  GSL_MULTIROOT_FN_EVAL_F_DF (FDF, x, f, J);

  for (i = 0; i < n; i++)
    {
      gsl_vector_set (dx, i, 0.0);
    }

  /* store column norms in diag */

  for (j = 0; j < n; j++)
    {
      double sum = 0;
      for (i = 0; i < n ; i++)
        {
          double Jij = gsl_matrix_get(J, i, j);
          sum += Jij * Jij ;
        }
      gsl_vector_set (diag, j, sqrt(sum));
    }

  /* set delta to factor |D x| or to factor if |D x| is zero */

  Dx = scaled_enorm (diag, x);

  factor = 100;

  state->delta = (Dx > 0) ? factor * Dx : factor ;

  /* Factorize J into QR decomposition */

  gsl_la_decomp_QR_impl (J, rdiag);

  gsl_la_unpack_QR_impl (J, rdiag, q, r);

  /* compute qtf = Q^T f */

  gsl_la_QTvec_QR_impl (J, f, qtf);

  return GSL_SUCCESS;
}

int
hybridj_iterate (void * vstate, gsl_multiroot_function_fdf * fdf, gsl_vector * x, gsl_vector * f, gsl_matrix * J, gsl_vector * dx)
{
  hybridj_state_t * state = (hybridj_state_t *) vstate;
  
  int signum ;

  size_t i;

  size_t n = fdf->n ;

  gsl_matrix * q = state->q;
  gsl_matrix * r = state->r;
  gsl_vector * rdiag = state->rdiag;
  gsl_vector * diag = state->diag;
  gsl_vector * qtf = state->qtf;

  /* compute dogleg step */

  dogleg (r, qtf, diag, delta, newton, gradient, dx);

  /* take a trial step */

  for (i = 0; i < n; i++)
    {
      double pi = gsl_vector_get (dx, i);
      double xi = gsl_vector_get (x, i);
      gsl_vector_set (x_trial, i, xi + pi);
    }

  pnorm = scaled_enorm (diag, dx);

  /* evaluate function at x + p */

  GSL_MULTIROOT_FN_EVAL_F (fdf, x, f);

  fnorm1 = enorm(f);

  
  /* compute the scaled actual reduction */


  if (fnorm1 < fnorm)
    {
      double u = fnorm1 / fnorm;
      actred = 1 - u * u;
    }
  else
    {
      actred = -1;
    }

  /* compute the scaled predicted reduction */
  
  fnorm1p = predicted_value (qtf, r, dx);

  if (fnorm1p < fnorm)
    {
      double u = fnorm1p / fnorm ;
      prered = 1 -  u * u ;
    }
  else
    {
      preref = 0;
    }


  if (prered > 0)
    {
      ratio = actred / prered ;
    }
  else
    {
      ratio = 0;
    }

  if (ratio < p1)
    {
      nsuc = 0;
      ncfail++;
      delta *= p5;
    }
  else
    {
      ncfail = 0;
      ncsuc++;
      if (ratio >= p5 || ncsuc > 1)
        delta = GSL_MAX(delta, pnorm/p5);
      if (fabs(ratio-1) <= p1) 
        delta = pnorm/p5;
    }

  if (ratio >= p0001)
    {
      /* successful iteration */


    }

  /* determine the progress of the iteration */
  

  


  /* rank-1 update of the jacobian */

  gsl_la_decomp_QR_impl (J, state->rdiag);

  gsl_la_unpack_QR_impl (J, state->rdiag, state->q, state->r);

  /* compute qtf = Q^T f */

  gsl_la_QTvec_QR_impl (J, f, state->qtf);

  return GSL_SUCCESS;
}


void
hybridj_free (void * vstate)
{
  hybridj_state_t * state = (hybridj_state_t *) vstate;
 
  gsl_vector_free(state->gradient);
  gsl_vector_free(state->newton);
  gsl_vector_free(state->qtf);
  gsl_vector_free(state->diag);
  gsl_vector_free(state->rdiag);
  gsl_matrix_free(state->r);
  gsl_matrix_free(state->q);
}

static const gsl_multiroot_fdfsolver_type hybridj_type =
{"hybridj",				/* name */
 sizeof (hybridj_state_t),
 &hybridj_alloc,
 &hybridj_set,
 &hybridj_iterate,
 &hybridj_free};

const gsl_multiroot_fdfsolver_type  * gsl_multiroot_fdfsolver_hybridj = &hybridj_type;
