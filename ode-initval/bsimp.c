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

  gsl_matrix      *  d;      /* workspace for extrapolation         */
  gsl_matrix      *  a_mat;  /* workspace for linear system matrix  */
  gsl_vector_int  *  p_vec;  /* workspace for LU permutation vector */

  double x[KMAXX];  /* workspace for extrapolation */

  /* state info */
  size_t kchoice;
  double hnext;
  double eps;

  /* workspace for extrapolation step */
  double * ysav;       
  double * yseq;
  double * yp;
  double * extr_work;
  double * dfdt;
  gsl_matrix * dfdy;

  /* workspace for the basic stepper */
  gsl_vector * rhs_temp_vec;
  gsl_vector * delta_temp_vec;
  double * delta;
};


/* Calculate a choice for the "order" of
 * the method, using the Deuflhard criteria.
 */
static size_t
bsimp_deuf_kchoice(double eps, size_t dimension)
{
  const double eps1 = SAFE1 * eps;

  double ab[KMAXX+1];
  double alf[KMAXX][KMAXX];

  size_t k, iq;

  ab[0] = nseq[0] + 1.0;

  for(k=0; k<KMAXX; k++) {
    ab[k + 1] = ab[k] + nseq[k + 1];
  }

  for(iq=1; iq<KMAXX; iq++) {
    for(k=0; k < iq-1; k++) {
      const double tmp1 = ab[k + 1] - ab[iq + 1];
      const double tmp2 = (ab[iq + 1] - ab[0] + 1.0) * (2*k + 1);
      alf[k][iq] = pow(eps1, tmp1/tmp2);
    }
  }

  ab[0] += dimension;

  for(k=0; k<KMAXX; k++) {
    ab[k + 1] = ab[k] + nseq[k + 1];
  }

  for(k=1; k < KMAXX-1; k++) {
    if(ab[k + 1] > ab[k] * alf[k - 1][k]) break;
  }

  return k;
}


static int bsimp_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);


/* Perform basic allocation of temporaries
 * given the dimension of the system vector.
 */
static void
bsimp_alloc(gsl_odeiv_step_bsimp * self, size_t dim)
{
  self->d     = gsl_matrix_alloc(KMAXX, dim);
  self->a_mat = gsl_matrix_alloc(dim, dim);
  self->p_vec = gsl_vector_int_alloc(dim);

  self->ysav = (double *) malloc(dim * sizeof(double));
  self->yseq = (double *) malloc(dim * sizeof(double));
  self->dfdt = (double *) malloc(dim * sizeof(double));
  self->yp   = (double *) malloc(dim * sizeof(double));
  self->extr_work = (double *) malloc(dim * sizeof(double));

  self->dfdy = gsl_matrix_alloc(dim, dim);

  self->delta_temp_vec = gsl_vector_alloc(dim);
  self->rhs_temp_vec = gsl_vector_alloc(dim);

  self->delta = (double *) malloc(dim * sizeof(double));
}


/* Perform basic de-allocation of temporaries.
 */
static void
bsimp_dealloc(gsl_odeiv_step_bsimp * self)
{
  if(self->delta != 0) free(self->delta);

  if(self->rhs_temp_vec != 0) gsl_vector_free(self->rhs_temp_vec);
  if(self->delta_temp_vec != 0) gsl_vector_free(self->delta_temp_vec);

  if(self->dfdy != 0) gsl_matrix_free(self->dfdy);

  if(self->extr_work != 0) free(self->extr_work);
  if(self->yp != 0) free(self->yp);
  if(self->dfdt != 0) free(self->dfdt);
  if(self->yseq != 0) free(self->yseq);
  if(self->ysav != 0) free(self->ysav);

  if(self->p_vec != 0) gsl_vector_int_free(self->p_vec);
  if(self->a_mat != 0) gsl_matrix_free(self->a_mat);
  if(self->d != 0) gsl_matrix_free(self->d);

  self->delta = 0;

  self->rhs_temp_vec = 0;
  self->delta_temp_vec = 0;

  self->dfdy = 0;

  self->extr_work = 0;
  self->yp = 0;
  self->dfdt = 0;
  self->yseq = 0;
  self->ysav = 0;

  self->p_vec = 0;
  self->a_mat = 0;
  self->d = 0;
}


static int
bsimp_reset(void * self)
{
  if(self != 0) {
    gsl_odeiv_step_bsimp * my = (gsl_odeiv_step_bsimp *) self;
    bsimp_dealloc(my);
    if(my->parent.dimension > 0) {
      bsimp_alloc(my, my->parent.dimension);
      my->kchoice = bsimp_deuf_kchoice(my->eps, my->parent.dimension);
      my->parent.order = 2 * my->kchoice;
    }
  }

  return GSL_SUCCESS;
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
gsl_odeiv_step_bsimp_new(double eps)
{  
  gsl_odeiv_step_bsimp * s = (gsl_odeiv_step_bsimp *) malloc(sizeof(gsl_odeiv_step_bsimp));
  if(s != 0) {
    gsl_odeiv_step_construct(&(s->parent),
      "bsimp",
      bsimp_step,
      bsimp_reset,
      bsimp_free,
      1,
      0,
      0);
      s->a_mat = 0;
      s->p_vec = 0;
      s->d = 0;
      s->extr_work = 0;
      s->yp = 0;
      s->ysav = 0;
      s->yseq = 0;
      s->dfdt = 0;
      s->dfdy = 0;
      s->eps = eps;
      s->rhs_temp_vec = 0;
      s->delta_temp_vec = 0;
      s->delta = 0;
  }
  return (gsl_odeiv_step *) s;
}



/*
   Use polynomial extrapolation to evaluate dim functions at x = 0 by fitting a 
   polynomial to a
   sequence of estimates with progressively smaller values x = xest, and 
   corresponding function
   vectors yest[1..dim]. This call is number iest in the sequence of calls.
   Extrapolated function
   values are output as y_0[1..dim], and their estimated error is output as y_0_err[1..dim].
 */
static void
pzextr(
  gsl_matrix * d,
  const double x[],
  int iest,
  double xest,
  double yest[],
  double y_0[],
  double y_0_err[],
  double work[],
  int dim)
{
  size_t j, k;

  memcpy(y_0_err, yest, dim * sizeof(double));
  memcpy(y_0,     yest, dim * sizeof(double));

  if(iest == 0) {
    for(j=0; j<dim; j++) {
      gsl_matrix_set(d, 0, j, yest[j]);
    }
  }
  else {
    memcpy(work, yest, dim * sizeof(double));

    for(k=0; k<iest; k++) {
      double delta = 1.0/(x[iest-k-1] - xest);
      const double f1 = delta * xest;
      const double f2 = delta * x[iest-k-1];

      /* Propagate another diagonal. */
      for(j=0; j<dim; j++) {
        const double q_kj = gsl_matrix_get(d, k, j);
        gsl_matrix_set(d, k, j, y_0_err[j]);
        delta = work[j] - q_kj;
        y_0_err[j] = f1 * delta;
        work[j]    = f2 * delta;
        y_0[j]    += y_0_err[j];
      }
    }

    for(j=0; j<dim; j++) {
      gsl_matrix_set(d, iest, j, y_0_err[j]);
    }
  }
}


#if 0
/*
   Exact substitute for pzextr, but uses diagonal rational function extrapolation instead of poly-
   nomial extrapolation.
 */
static void 
rzextr(
  gsl_matrix * d,
  const double x[],
  int iest,
  double xest,
  double yest[],
  double y_0[],
  double y_0_err[],
  double work[],
  int dim)
{
  size_t k, j;

  if(iest == 1) {
    for(j=0; j<dim; j++) {
      y_0[j] = yest[j];
      y_0_err[j] = yest[j];
      gsl_matrix_set(d, 0, j, yest[j]);
    }
  }
  else {
    for(k=0; k<iest; k++) {
      work[k + 1] = x[iest - k] / xest;
    }

    /* Next diagonal. */
    for(j=0; j<dim; j++) {
      double ddy;
      double v  = gsl_matrix_get(d, 0, j);
      double yy = yest[j];
      double c  = yest[j];
      gsl_matrix_set(d, 0, j, yest[j]);

      for(k=1; k<=iest; k++) {
        double b1 = work[k] * v;
        double b = b1 - c;
        if(b != 0) {
          b = (c - v) / b;
          ddy = c * b;
          c = b1 * b;
        }
        else {
          /* Care needed to avoid division by 0. */
          ddy = v;
	}
        if(k != iest) {
          v = gsl_matrix_get(d, k, j);
	}
	gsl_matrix_set(d, k, j, ddy);
        yy += ddy;
      }
      y_0_err[j] = ddy;
      y_0[j] = yy;
    }
  }
}
#endif /* 0 */


/* Basic implicit Bulirsch-Stoer step.
 * Divide the step h_total into nstep smaller
 * steps and do the Bader-Dueflhard
 * semi-implicit iteration.
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
  const double h = h_total / nstep;

  double t = t0 + h;

  int signum;
  unsigned int i, j;
  unsigned int n_inter;

  /* shameless use of somebody else's temp work space */
  double * ytemp = step->extr_work;
  gsl_vector ytemp_vec;
  ytemp_vec.data = ytemp;
  ytemp_vec.size = step->parent.dimension;
  ytemp_vec.parent = 0;
  ytemp_vec.stride = 1;

  /* Calculate the matrix for the linear system. */
  for(i=0; i<dim; i++) {
    for(j=0; j<dim; j++) {
      gsl_matrix_set(step->a_mat, i, j, -h * gsl_matrix_get(dfdy, i, j));
    }
    gsl_matrix_set(step->a_mat, i, i, gsl_matrix_get(step->a_mat, i, i) + 1.0);
  }

  /* LU decomposition for the linear system. */
  gsl_la_decomp_LU_impl(step->a_mat, step->p_vec, &signum);


  /* Initial step. */

  for(i=0; i<dim; i++) {
    ytemp[i] = h * (yp[i] + h * dfdt[i]);
  }

  gsl_la_solve_LU_impl(step->a_mat, step->p_vec, &ytemp_vec, step->delta_temp_vec);

  for (i=0; i<dim; i++) {
    const double di = gsl_vector_get(step->delta_temp_vec, i);
    step->delta[i] = di;
    ytemp[i] = y[i] + di;
  }


  /* Intermediate steps. */

  GSL_ODEIV_FN_EVAL(sys, t, ytemp, yout);

  for (n_inter=1; n_inter<nstep; n_inter++) {

    for (i=0; i<dim; i++) {
      gsl_vector_set(step->rhs_temp_vec, i, h * yout[i] - step->delta[i]);
    }

    gsl_la_solve_LU_impl(step->a_mat, step->p_vec, step->rhs_temp_vec, step->delta_temp_vec);

    for(i=0; i<dim; i++) {
      step->delta[i] += 2.0 * gsl_vector_get(step->delta_temp_vec, i);
      ytemp[i] += step->delta[i];
    }

    t += h;

    GSL_ODEIV_FN_EVAL(sys, t, ytemp, yout);
  }


  /* Final step. */

  for(i=0; i<dim; i++) {
    gsl_vector_set(step->rhs_temp_vec, i, h * yout[i] - step->delta[i]);
  }

  gsl_la_solve_LU_impl(step->a_mat, step->p_vec, step->rhs_temp_vec, step->delta_temp_vec);

  for(i=0; i<dim; i++) {
    yout[i] = ytemp[i] + gsl_vector_get(step->delta_temp_vec, i);
  }

  return GSL_SUCCESS;
}


/* Perform the basic semi-implicit extrapolation
 * step at a Deuflhard determined order.
 */
static int
bsimp_step(
  void * self,
  double t,
  double h,
  double y[],
  double yerr[],
  const double dydt_in[],
  double dydt_out[],
  const gsl_odeiv_system * sys)
{
  gsl_odeiv_step_bsimp * my = (gsl_odeiv_step_bsimp *) self;

  const double t_local = t;

  size_t k;

  if(sys->jacobian == 0) {
    return GSL_EFAULT; /* FIXME: error condition */
  }

  if(sys->dimension <= 0) {
    return GSL_EINVAL;
  }

  if(h + t_local == t_local) {
    return GSL_EUNDRFLW; /* FIXME: error condition */
  }

  /* Trap and perform state initialization,
   * which allocates work space and calculates
   * the effective "order" of the method.
   */
  if(my->parent.dimension != sys->dimension) {
    my->parent.dimension = sys->dimension;
    bsimp_dealloc(my);
    bsimp_alloc(my, my->parent.dimension);
    my->kchoice = bsimp_deuf_kchoice(my->eps, my->parent.dimension);
    my->parent.order = 2 * my->kchoice;
    my->hnext   = -GSL_SQRT_DBL_MAX;
  }

  memcpy(my->ysav, y, my->parent.dimension * sizeof(double));

  /* Evaluate the derivative. */
  if(dydt_in != 0) {
    memcpy(my->yp, dydt_in, my->parent.dimension * sizeof(double));
  }
  else {
    GSL_ODEIV_FN_EVAL(sys, t_local, y, my->yp);
  }

  /* Evaluate the Jacobian for the system. */
  GSL_ODEIV_JA_EVAL(sys, t_local, y, my->dfdy->data, my->dfdt);

  /* Make a series of refined extrapolations,
   * up the specified maximum order, which
   * was calculated based on the Deuflhard
   * criterion upon state initialization.
   */
  for(k=0; k <= my->kchoice; k++) {

    const double xest = (h/nseq[k]) * (h/nseq[k]);

    bsimp_step_local(self,
  		     my->ysav,
		     my->yp,
  		     my->dfdt,
  		     my->dfdy,
  		     my->parent.dimension,
  		     t_local,
  		     h,
  		     nseq[k],
  		     my->yseq,
  		     sys);

    my->x[k] = xest;
    pzextr(my->d, my->x, k, xest, my->yseq, y, yerr, my->extr_work, my->parent.dimension);
  }

  return GSL_SUCCESS;
}
