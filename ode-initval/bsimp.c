/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_linalg.h>
#include "odeiv_util.h"
#include "gsl_odeiv.h"

#define KMAXX 7
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1


static const int nseq[IMAXX + 1] = { 0, 2, 6, 10, 14, 22, 34, 50, 70 };


typedef  struct gsl_odeiv_step_bsimp_struct  gsl_odeiv_step_bsimp;

struct gsl_odeiv_step_bsimp_struct
{
  gsl_odeiv_step parent;  /* inherits from gsl_odeiv_step */

  gsl_matrix  *  a;
  gsl_vector  *  indx;

  double ** d;
  double x[KMAXX];

  double ab[IMAXX + 1];
  double alf[KMAXX + 1][KMAXX + 1];

  double err[KMAXX];

  int first;
  int kmax;
  int kopt;
  int nvold;

  double epsold;
  double xnew;

  double * ysav;
  double * yseq;

  double  * dfdx;
  double ** dfdy;

  gsl_vector rhs_temp_vec;
  gsl_vector delta_temp_vec;
};


static int  bsimp_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
static int  bsimp_reset(void * self);
static void bsimp_free(void * self);


gsl_odeiv_step *
gsl_odeiv_step_bsimp_new()
{
  
  gsl_odeiv_step_bsimp * s = (gsl_odeiv_step_bsimp *) malloc(sizeof(gsl_odeiv_step_bsimp));
  if(s != 0) {
    gsl_odeiv_step_construct(&(s->parent),
      "bsimp",
      bsimp_step,
      bsimp_reset,
      bsimp_free,
      0,
      0,
      2);
      s->a = 0;
      s->indx = 0;
      s->d = 0;
      s->first = 1;
      s->nvold = -1;
      s->epsold = -1.0;
      s->ysav = 0;
      s->yseq = 0;
      s->dfdx = 0;
      s->dfdy = 0;
  }
  return (gsl_odeiv_step *) s;
}




extern double **d, *x;		/* Defined in bsstep. */

void 
pzextr (int iest, double xest, double yest[], double yz[], double dy[], int nv)
/*
   Use polynomial extrapolation to evaluate nv functions at x = 0 by tting a polynomial to a
   sequence of estimates with progressively smaller values x = xest, and corresponding function
   vectors yest[1..nv]. This call is number iest in the sequence of calls. Extrapolated function
   values are output as yz[1..nv], and their estimated error is output as dy[1..nv].
 */
{
  int k1, j;
  double q, f2, f1, delta, *c;
  c = vector (1, nv);
  x[iest] = xest;		/* Save current independent variable. */
  for (j = 1; j <= nv; j++)
    dy[j] = yz[j] = yest[j];
  if (iest == 1)
    {				/* Store first estimate in first column. */
      for (j = 1; j <= nv; j++)
	d[j][1] = yest[j];
    }
  else
    {
      for (j = 1; j <= nv; j++)
	c[j] = yest[j];
      for (k1 = 1; k1 < iest; k1++)
	{
	  delta = 1.0 / (x[iest - k1] - xest);
	  f1 = xest * delta;
	  f2 = x[iest - k1] * delta;
	  for (j = 1; j <= nv; j++)
	    {			/* Propagate tableau 1 diagonal more. */
	      q = d[j][k1];
	      d[j][k1] = dy[j];
	      delta = c[j] - q;
	      dy[j] = f1 * delta;
	      c[j] = f2 * delta;
	      yz[j] += dy[j];
	    }
	}
      for (j = 1; j <= nv; j++)
	d[j][iest] = dy[j];
    }
  free_vector (c, 1, nv);
}



extern double **d, *x; /* Defined in bsstep. */

void 
rzextr (int iest, double xest, double yest[], double yz[], double dy[], int nv)
/*
   Exact substitute for pzextr, but uses diagonal rational function extrapolation instead of poly-
   nomial extrapolation.
 */
{
  int k, j;
  double yy, v, ddy, c, b1, b, *fx;
  fx = vector (1, iest);
  x[iest] = xest;		/* Save current independent variable. */
  if (iest == 1)
    for (j = 1; j <= nv; j++)
      {
	yz[j] = yest[j];
	d[j][1] = yest[j];
	dy[j] = yest[j];
      }
  else
    {
      for (k = 1; k < iest; k++)
	fx[k + 1] = x[iest - k] / xest;
      for (j = 1; j <= nv; j++)
	{			/* Evaluate next diagonal in tableau. */
	  v = d[j][1];
	  d[j][1] = yy = c = yest[j];
	  for (k = 2; k <= iest; k++)
	    {
	      b1 = fx[k] * v;
	      b = b1 - c;
	      if (b)
		{
		  b = (c - v) / b;
		  ddy = c * b;
		  c = b1 * b;
		}
	      else		/* Care needed to avoid division by 0. */
		ddy = v;
	      if (k != iest)
		v = d[j][k];
	      d[j][k] = ddy;
	      yy += ddy;
	    }
	  dy[j] = ddy;
	  yz[j] = yy;
	}
    }
  free_vector (fx, 1, iest);
}


static int
bsimp_step_local(
  gsl_odeiv_step_bsimp * step,
  const double y[],
  const double yp[],
  const double dfdt[],
  const gsl_matrix * dfdy,
  unsigned int dim,
  double t0,
  double h_total,
  unsigned int nstep,
  double yout[],
  const gsl_odeiv_system * sys
)
/* simpr (double y[], double dydx[], double dfdx[], double **dfdy, int n,
       double xs, double htot, int nstep, double yout[],
       void (*derivs) (double, double[], double[]))
       */
{
  int signum;
  unsigned int i, j;
  unsigned int n_inter;

  const double h = h_total / nstep;

  double t = t0 + h;

  /* shameless confusion... use these for temp space */
  double * ytemp = step->ysav;
  double * delta_temp = step->yseq;
  double * delta = step->dfdx;

  step->delta_temp_vec.data = delta_temp;
  step->delta_temp_vec.size = dim;
  step->rhs_temp_vec.data = 0 /* FIXME */;
  step->rhs_temp_vec.size = dim;

  for(i=0; i<dim; i++) {
    for(j=0; j<dim; j++) {
      gsl_matrix_set(step->a, i, j, -h * gsl_matrix_get(dfdy, i, j));
    }
    gsl_matrix_set(step->a, i, i, gsl_matrix_get(step->a, i, i) + 1.0);
  }

  /* LU decomposition for the linear system. */
  gsl_la_decomp_LU_impl(step->a, step->indx, &signum);


  /* Initial step. */

  for(i=0; i<dim; i++) {
    ytemp[i] = h * (yp[i] + h * dfdt[i]);
  }

  gsl_la_solve_LU_impl(step->a, step->indx, ytemp, &step->delta_temp_vec);

  for (i=0; i<dim; i++) {
    delta[i] = delta_temp[i];
    ytemp[i] = y[i] + delta_temp[i];
  }


  /* Intermediate steps. */

  GSL_ODEIV_FN_EVAL(sys, t, ytemp, yout);

  for (n_inter=1; n_inter<nstep; n_inter++) {

    for (i=0; i<dim; i++) {
      step->rhs_temp_vec.data[i] = h * yout[i] - delta[i];
    }

    gsl_la_solve_LU_impl(step->a, step->indx, &step->rhs_temp_vec, &step->delta_temp_vec);

    for(i=0; i<dim; i++) {
      delta[i] += 2.0 * delta_temp[i];
      ytemp[i] += delta[i];
    }

    t += h;

    GSL_ODEIV_FN_EVAL(sys, t, ytemp, yout);
  }


  /* Final step. */

  for(i=0; i<dim; i++) {
    step->rhs_temp_vec.data[i] = h * yout[i] - delta[i];
  }

  gsl_la_solve_LU_impl(step->a, step->indx, &step->rhs_temp_vec, &step->delta_temp_vec);

  for(i=0; i<dim; i++) {
    yout[i] = ytemp[i] + delta_temp[i];
  }

  return GSL_SUCCESS;
}




static int
bsimp_step(
  void * self,
  double t,
  double htry,
  double y[],
  double yerr[],
  const double dydt_in[],
  double dydt_out[],
  const gsl_odeiv_system * sys)

/*
void
stifbs (double y[], double dydx[], int nv, double *xx, double htry,
        double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs) (double, double[], double[]))
	*/
/*
   Semi - implicit extrapolation step for integrating stiff o.d.e.'s, with 
   monitoring of local truncation
   error to adjust stepsize.Input are the dependent variable vector y[1..n]
   and its derivative
   dydx[1..n] at the starting value of the independent variable x.Also input
   are the stepsize to
   be attempted htry, the required accuracy eps, and the vector yscal[1..n]
   against which
   the error is scaled.On output, y and x are replaced by their new values,
   hdid is the stepsize
   that was actually accomplished, and hnext is the estimated next
   stepsize.
   derivs is a user - nsupplied routine that computes the derivatives
   of the right - hand side with respect to x, while
   jacobn (a          xed name)
   is a user - supplied routine that computes the Jacobi matrix of
   derivatives
   of the right - hand side with respect to the components of y.Be
   sure to set htry on successive
   steps to the value of hnext returned from the previous step,
   as is the case if the routine is
   called by odeint.
 */
{
  gsl_odeiv_step_bsimp * my = (gsl_odeiv_step_bsimp *) self;

  size_t dim = sys->dimension;

  int i, iq, k, kk, km;

  double eps1, errmax, fact, red, scale, work, wrkmin, xest;


  int reduct;

  int exitflag = 0;

  double h = htry;


  if(sys->dimension <= 0) {
    return GSL_EINVAL;
  }


/*
  d = matrix (1, nv, 1, KMAXX);

*/
/*
  dfdx = vector (1, nv);
  dfdy = matrix (1, nv, 1, nv);
*/
  /*
  err = vector (1, KMAXX);
  x = vector (1, KMAXX);
  */
  /* yerr = vector (1, nv); */

/*
  ysav = vector (1, nv);
  yseq = vector (1, nv);
*/

  if (eps != my->epsold || dim != my->nvold)
    {
      /* Reinitialize also if nv has changed. */

      *hnext = my->xnew = -1.0e29;
      eps1 = SAFE1 * eps;
      my->ab[1] = nseq[1] + 1;

      for (k = 1; k <= KMAXX; k++)
	my->ab[k + 1] = my->ab[k] + nseq[k + 1];
 
      for (iq = 2; iq <= KMAXX; iq++)
	{
	  for (k = 1; k < iq; k++)
	    my->alf[k][iq] = pow (eps1, ((my->ab[k + 1] - my->ab[iq + 1]) /
				 ((my->ab[iq + 1] - my->ab[1] + 1.0) * (2 * k + 1))));
	}
      my->epsold = eps;
      my->nvold = dim;  /* Save nv. */
      my->ab[1] += dim; /* Add cost of Jacobian evaluations to work coefficients. */
      for (k = 1; k <= KMAXX; k++)
	my->ab[k + 1] = my->ab[k] + nseq[k + 1];

      /* determine an optimal order */
      for (my->kopt = 2; my->kopt < KMAXX; my->kopt++)
	if (my->ab[my->kopt + 1] > my->ab[my->kopt] * my->alf[my->kopt - 1][my->kopt])
	  break;
      my->kmax = my->kopt;
    }


  for(i=0; i<dim; i++)
    my->ysav[i] = y[i];

/*
  jacobn (*xx, y, dfdx, dfdy, dim);
*/
  if (*xx != my->xnew || h != (*hnext))
    {
      my->first = 1;
      my->kopt = my->kmax;
    }

  reduct = 0;

  while(1) {

    /* tries extrapolations up to a max order... */

    for (k = 1; k <= my->kmax; k++)
      {
        my->xnew = (*xx) + h;

        if (my->xnew == (*xx))
          nrerror ("step size underflow in stifbs");

/*
        simpr (my->ysav, dydx, dfdx, dfdy, nv, *xx, h, nseq[k], my->yseq, derivs);
*/
        xest = SQR (h / nseq[k]);


        /* The rest of the routine is identical to
           bsstep. */

        pzextr (k, xest, my->yseq, y, yerr, dim);


        if (k != 1)
          {
            errmax = TINY;
            for(i=0; i<dim; i++) {
              errmax = FMAX (errmax, fabs (yerr[i] / yscal[i]));
            }
            errmax /= eps;
            km = k - 1;
            my->err[km] = pow (errmax / SAFE1, 1.0 / (2 * km + 1));
          }

        if (k != 1 && (k >= my->kopt - 1 || my->first))
          {
            if (errmax < 1.0)
              {
        	exitflag = 1;
        	break;
              }
            if (k == my->kmax || k == my->kopt + 1)
              {
        	red = SAFE2 / my->err[km];
        	break;
              }
            else if (k == my->kopt && my->alf[my->kopt - 1][my->kopt] < my->err[km])
              {
        	red = 1.0 / my->err[km];
        	break;
              }
            else if (my->kopt == my->kmax && my->alf[km][my->kmax - 1] < my->err[km])
              {
        	red = my->alf[km][my->kmax - 1] * SAFE2 / my->err[km];
        	break;
              }
            else if (my->alf[km][my->kopt] < my->err[km])
              {
        	red = my->alf[km][my->kopt - 1] / my->err[km];
        	break;
              }
          }
      }

    if (exitflag) break;

    red = FMIN (red, REDMIN);
    red = FMAX (red, REDMAX);
    h *= red;
    reduct = 1;
  }

  *xx = my->xnew;
  *hdid = h;
  my->first = 0;
  wrkmin = 1.0e35;
  for (kk = 1; kk <= km; kk++)
    {
      fact = FMAX (my->err[kk], SCALMX);
      work = fact * my->ab[kk + 1];
      if (work < wrkmin)
	{
	  scale = fact;
	  wrkmin = work;
	  my->kopt = kk + 1;
	}
    }

  *hnext = h / scale;

  if (my->kopt >= k && my->kopt != my->kmax && !reduct)
    {
      fact = FMAX (scale / my->alf[my->kopt - 1][my->kopt], SCALMX);
      if (my->ab[my->kopt + 1] * fact <= wrkmin)
	{
	  *hnext = h / fact;
	  my->kopt++;
	}
    }

/*
  free_vector (yseq, 1, nv);
  free_vector (ysav, 1, nv);
  free_vector (yerr, 1, nv);
  free_vector (x, 1, KMAXX);
  free_vector (err, 1, KMAXX);
  free_matrix (dfdy, 1, nv, 1, nv);
  free_vector (dfdx, 1, nv);
  free_matrix (d, 1, nv, 1, KMAXX);
  */
}
