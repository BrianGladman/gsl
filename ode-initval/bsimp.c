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


static const int nseq[IMAXX] = { 2, 6, 10, 14, 22, 34, 50, 70 };


typedef  struct gsl_odeiv_step_bsimp_struct  gsl_odeiv_step_bsimp;

struct gsl_odeiv_step_bsimp_struct
{
  gsl_odeiv_step parent;  /* inherits from gsl_odeiv_step */

  size_t dimension;
  size_t dimension_old;

  gsl_matrix  *  a;
  gsl_vector  *  indx;

  gsl_matrix  * d;
  double x[KMAXX];
  double err[KMAXX];

  double ab[IMAXX];
  double alf[KMAXX][KMAXX];

  int first;
  int kmax;
  int kopt;

  double xnew;

  double hdid;
  double hnext;

  const double * yscal;
  double eps;
  double epsold;

  double * ysav;
  double * yseq;

  double * dfdt;
  gsl_matrix * dfdy;

  gsl_vector rhs_temp_vec;
  gsl_vector delta_temp_vec;
};


static int  bsimp_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);


/* Perform basic allocation of temporaries
 * given the dimension of the system vector.
 */
static int
bsimp_alloc(gsl_odeiv_step_bsimp * self, size_t dim)
{
  self->ysav = (double *) malloc(dim * sizeof(double));
  self->yseq = (double *) malloc(dim * sizeof(double));
  self->dfdt = (double *) malloc(dim * sizeof(double));

  self->dfdy = gsl_matrix_alloc(dim, dim);
  self->d = gsl_matrix_double_alloc(dim, KMAXX);
  self->a = gsl_matrix_double_alloc(dim, dim);
}


/* Perform basic dallocation of temporaries.
 */
static void
bsimp_dealloc(gsl_odeiv_step_bsimp * self)
{
  if(self->ysav != 0) free(self->ysav);
  if(self->yseq != 0) free(self->yseq);
  if(self->dfdt != 0) free(self->dfdt);
  if(self->dfdy != 0) gsl_matrix_free(self->dfdy);
  if(self->d != 0) gsl_matrix_free(self->d);
  if(self->a != 0) gsl_matrix_free(self->a);
}


static void
bsimp_reset(void * self)
{
  if(self != 0) {
    gsl_odeiv_step_bsimp * my = (gsl_odeiv_step_bsimp *) self;
    my->first = 1;
    my->dimension_old = -1;
    my->epsold = -1.0;
    bsimp_dealloc(my);
    if(my->dimension > 0) {
      bsimp_alloc(my, my->dimension);
    }
  }
}


static void
bsimp_free(void * self)
{
  if(self != 0) {
    bsimp_dealloc((gsl_odeiv_step_bsimp *) self);
    free(self);
  }
}


gsl_odeiv_step *
gsl_odeiv_step_bsimp_new(const double * yscal, double eps)
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
      2 /* FIXME */);
      s->a = 0;
      s->indx = 0;
      s->d = 0;
      s->ysav = 0;
      s->yseq = 0;
      s->dfdt = 0;
      s->dfdy = 0;
      s->yscal = yscal; /* FIXME */
      s->eps = eps;
      s->dimension = 0;
  }
  return (gsl_odeiv_step *) s;
}


/* Reset the Deuflhard error control data.
 */
static void
bsimp_reset_deuf(gsl_odeiv_step_bsimp * self)
{
  size_t k;
  size_t iq;

  /* my->hnext = my->xnew = -1.0e29; */
  const double eps1 = SAFE1 * self->eps;

  self->ab[0] = nseq[0] + 1;

  for(k=0; k<KMAXX; k++) {
    self->ab[k + 1] = self->ab[k] + nseq[k + 1];
  }

  for(iq=1; iq < KMAXX; iq++) {
    for(k=0; k < iq; k++) {
      const double tmp1 = self->ab[k + 1] - self->ab[iq + 1];
      const double tmp2 = (self->ab[iq + 1] - self->ab[1] + 1.0) * (2 * k + 1);
      self->alf[k][iq] = pow(eps1, tmp1/tmp2);
    }
  }

  self->epsold = self->eps;
  self->ab[0] += self->dimension;
  self->dimension_old = self->dimension;

  for(k=0; k < KMAXX; k++) {
    self->ab[k + 1] = self->ab[k] + nseq[k + 1];
  }

  /* determine an optimal order */
  for(k = 1; k < KMAXX-1; k++) {
    if (self->ab[k + 1] > self->ab[k] * self->alf[k - 1][k]) {
      self->kopt = k;
      break;
    }
  }
  self->kmax = self->kopt;
}


void 
pzextr(double ** d, double * x, int iest, double xest, double yest[], double yz[], double dy[], int nv)
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



void 
rzextr(double ** d, double * x, int iest, double xest, double yest[], double yz[], double dy[], int nv)
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


/* Basic implicit Bulirsch-Stoer step.
 */
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
  const gsl_odeiv_system * sys)
{
  int signum;
  unsigned int i, j;
  unsigned int n_inter;

  const double h = h_total / nstep;

  double t = t0 + h;

  /* shameless confusion... use these for temp space */
  double * ytemp = step->ysav;
  double * delta_temp = step->yseq;
  double * delta = step->dfdt;

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
bsimp_iter_test(
  gsl_odeiv_step_bsimp * self,
  size_t k,
  double * reduction_factor,
  const double * yerr,
  const double * yscal)
{
  const size_t km = k - 1;
  double errmax = DBL_MIN;
  size_t i;

  for(i=0; i<self->dimension; i++) {
    errmax = GSL_MAX_DBL(errmax, fabs(yerr[i])/(fabs(yscal[i]) + GSL_SQRT_DBL_EPSILON));
  }
  errmax /= self->eps;

  self->err[km] = pow(errmax / SAFE1, 1.0 / (2 * km + 1));

  if(k >= self->kopt - 1 || self->first) {

    if(errmax < 1.0) {
      return 1; /* done */
    }
    else if (k == self->kmax || k == self->kopt + 1) {
      *reduction_factor = SAFE2 / self->err[km];
      return 0;
    }
    else if (k == self->kopt && self->alf[self->kopt - 1][self->kopt] < self->err[km]) {
      *reduction_factor = 1.0 / self->err[km];
      return 0;
    }
    else if (self->kopt == self->kmax && self->alf[km][self->kmax - 1] < self->err[km]) {
      *reduction_factor = self->alf[km][self->kmax - 1] * SAFE2 / self->err[km];
      return 0;
    }
    else if (self->alf[km][self->kopt] < self->err[km]) {
      *reduction_factor = self->alf[km][self->kopt - 1] / self->err[km];
      return 0;
    }
    else {
      return 0; /* FIXME: ??? */
    }

  }
  else {
    return 0; /* FIXME: ??? */
  }
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

  const double * yscal = ( my->yscal != 0 ? my->yscal : y );

  double h = htry;
  double t_local = t;

  int exitflag = 0;
  int k, k_exit;
  int reduct;

  double scale, wrkmin;

  if(sys->dimension <= 0) {
    return GSL_EINVAL;
  }

  if(my->dimension != sys->dimension) {
    my->dimension = sys->dimension;
    bsimp_dealloc(my);
    bsimp_alloc(my, my->dimension);
    my->hnext = my->xnew = -GSL_SQRT_DBL_MAX;
    bsimp_reset_deuf(my);
  }
  else if(my->eps != my->epsold) {
    my->hnext = my->xnew = -GSL_SQRT_DBL_MAX;
    bsimp_reset_deuf(my);
  }

  memcpy(my->ysav, y, my->dimension * sizeof(double));

  GSL_ODEIV_JA_EVAL(sys, t_local, y, my->dfdy, my->dfdt);

  if(t_local != my->xnew || h != my->hnext) {
    my->first = 1;
    my->kopt = my->kmax;
  }

  reduct = 0;

  while(!exitflag) {
    double xest;
    double red;

    /* tries extrapolations up to a max order... */

    for(k=1; k <= my->kmax && !exitflag; k++) {
      my->xnew = t_local + h;

      if (my->xnew == t_local)
        nrerror ("step size underflow in stifbs");

      bsimp_step_local(self,
        	       my->ysav, dydt_in,
        	       my->dfdt,
        	       my->dfdy,
        	       my->dimension,
        	       t_local,
        	       h,
        	       nseq[k],
        	       my->yseq,
        	       sys);

      xest = SQR (h / nseq[k]);

      /* The rest of the routine is identical to bsstep. */

      pzextr(my->d, my->x, k, xest, my->yseq, y, yerr, my->dimension);

      if(k != 1) {
        exitflag = bsimp_iter_test(my, k, &red, yerr, my->yscal);
      }

      k_exit = k;
    }

    if(!exitflag) {
      red = GSL_MIN_DBL(red, REDMIN);
      red = GSL_MAX_DBL(red, REDMAX);
      h *= red;
      reduct = 1;
    }
  }

  t_local = my->xnew;
  my->hdid = h;
  my->first = 0;
  wrkmin = GSL_DBL_MAX;

  for(k=0; k < k_exit-1; k++) {
    double fact = FMAX (my->err[k], SCALMX);
    double work = fact * my->ab[k + 1];
    if(work < wrkmin) {
        scale = fact;
        wrkmin = work;
        my->kopt = k + 1;
    }
  }

  my->hnext = h / scale;

  if(my->kopt >= k && my->kopt != my->kmax && !reduct) {
    double fact = FMAX (scale / my->alf[my->kopt - 1][my->kopt], SCALMX);
    if(my->ab[my->kopt + 1] * fact <= wrkmin) {
      my->hnext = h / fact;
      my->kopt++;
    }
  }

}
