/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Bader-Deuflhard implicit extrapolative stepper.
 * [Numer. Math., 41, 373 (1983)]
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include "odeiv_util.h"
#include "gsl_odeiv.h"


#define SEQUENCE_COUNT 8
#define SEQUENCE_MAX   7

/* Bader-Deuflhard extrapolation sequence */
static const int bd_sequence[SEQUENCE_COUNT] = { 2, 6, 10, 14, 22, 34, 50, 70 };



typedef  struct gsl_odeiv_step_bsimp_struct  gsl_odeiv_step_bsimp;

struct gsl_odeiv_step_bsimp_struct
{
  gsl_odeiv_step parent;     /* inherits from gsl_odeiv_step */

  gsl_matrix      *  d;      /* workspace for extrapolation         */
  gsl_matrix      *  a_mat;  /* workspace for linear system matrix  */
  gsl_permutation *  p_vec;  /* workspace for LU permutation        */

  double x[SEQUENCE_MAX];    /* workspace for extrapolation */

  /* state info */
  size_t k_choice;
  double h_next;
  double eps;

  /* workspace for extrapolation step */
  double * y_extrap_save;
  double * y_extrap_sequence;
  double * yp;
  double * extrap_work;
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
  const double safety_f  = 0.25;
  const double small_eps = safety_f * eps;

  double a_work[SEQUENCE_COUNT];
  double alpha[SEQUENCE_MAX][SEQUENCE_MAX];

  int i, k;

  a_work[0] = bd_sequence[0] + 1.0;

  for(k=0; k<SEQUENCE_MAX; k++) {
    a_work[k + 1] = a_work[k] + bd_sequence[k + 1];
  }

  for(i=0; i<SEQUENCE_MAX; i++) {
    alpha[i][i] = 1.0;
    for(k=0; k<i; k++) {
      const double tmp1 = a_work[k + 1] - a_work[i + 1];
      const double tmp2 = (a_work[i + 1] - a_work[0] + 1.0) * (2*k + 1);
      alpha[k][i] = pow(small_eps, tmp1/tmp2);
    }
  }

  a_work[0] += dimension;

  for(k=0; k<SEQUENCE_MAX; k++) {
    a_work[k + 1] = a_work[k] + bd_sequence[k + 1];
  }

  for(k=0; k < SEQUENCE_MAX-1; k++) {
    if(a_work[k + 2] > a_work[k+1] * alpha[k][k+1]) break;
  }

  return k;
}


static int bsimp_step(void * self, double t, double h, double y[], double y_err[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);


/* Perform basic allocation of temporaries
 * given the dimension of the system vector.
 */
static void
bsimp_alloc(gsl_odeiv_step_bsimp * self, size_t dim)
{
  self->d     = gsl_matrix_alloc(SEQUENCE_MAX, dim);
  self->a_mat = gsl_matrix_alloc(dim, dim);
  self->p_vec = gsl_permutation_alloc(dim);

  self->y_extrap_save = (double *) malloc(dim * sizeof(double));
  self->y_extrap_sequence = (double *) malloc(dim * sizeof(double));
  self->dfdt = (double *) malloc(dim * sizeof(double));
  self->yp   = (double *) malloc(dim * sizeof(double));
  self->extrap_work = (double *) malloc(dim * sizeof(double));

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

  if(self->extrap_work != 0) free(self->extrap_work);
  if(self->yp != 0) free(self->yp);
  if(self->dfdt != 0) free(self->dfdt);
  if(self->y_extrap_sequence != 0) free(self->y_extrap_sequence);
  if(self->y_extrap_save != 0) free(self->y_extrap_save);

  if(self->p_vec != 0) gsl_permutation_free(self->p_vec);
  if(self->a_mat != 0) gsl_matrix_free(self->a_mat);
  if(self->d != 0) gsl_matrix_free(self->d);

  self->delta = 0;

  self->rhs_temp_vec = 0;
  self->delta_temp_vec = 0;

  self->dfdy = 0;

  self->extrap_work = 0;
  self->yp = 0;
  self->dfdt = 0;
  self->y_extrap_sequence = 0;
  self->y_extrap_save = 0;

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
      my->k_choice = bsimp_deuf_kchoice(my->eps, my->parent.dimension);
      my->parent.order = 2 * my->k_choice;
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
      s->extrap_work = 0;
      s->yp = 0;
      s->y_extrap_save = 0;
      s->y_extrap_sequence = 0;
      s->dfdt = 0;
      s->dfdy = 0;
      s->eps = eps;
      s->rhs_temp_vec = 0;
      s->delta_temp_vec = 0;
      s->delta = 0;
  }
  return (gsl_odeiv_step *) s;
}



static void
poly_extrap(
  gsl_matrix * d,
  const double x[],
  const int i_step,
  const double x_i,
  const double y_i[],
  double y_0[],
  double y_0_err[],
  double work[],
  const size_t dim)
{
  size_t j, k;

  memcpy(y_0_err, y_i, dim * sizeof(double));
  memcpy(y_0,     y_i, dim * sizeof(double));

  if(i_step == 0) {
    for(j=0; j<dim; j++) {
      gsl_matrix_set(d, 0, j, y_i[j]);
    }
  }
  else {
    memcpy(work, y_i, dim * sizeof(double));

    for(k=0; k<i_step; k++) {
      double delta = 1.0/(x[i_step-k-1] - x_i);
      const double f1 = delta * x_i;
      const double f2 = delta * x[i_step-k-1];

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
      gsl_matrix_set(d, i_step, j, y_0_err[j]);
    }
  }
}


/* Basic implicit Bulirsch-Stoer step.
 * Divide the step h_total into n_step smaller
 * steps and do the Bader-Deuflhard
 * semi-implicit iteration.
 */
static int
bsimp_step_local(
  gsl_odeiv_step_bsimp * step,
  const double y[],
  const double yp[],
  const double dfdt[],
  const gsl_matrix * dfdy,
  const size_t dim,
  const double t0,
  const double h_total,
  const unsigned int n_step,
  double y_out[],
  const gsl_odeiv_system * sys)
{
  const double h = h_total / n_step;

  double t = t0 + h;

  int signum;
  size_t i, j;
  size_t n_inter;

  /* shameless use of somebody else's temp work space */
  double * ytemp = step->extrap_work;
  gsl_vector ytemp_vec;
  ytemp_vec.data = ytemp;
  ytemp_vec.size = step->parent.dimension;
  ytemp_vec.block = 0;
  ytemp_vec.stride = 1;

  /* Calculate the matrix for the linear system. */
  for(i=0; i<dim; i++) {
    for(j=0; j<dim; j++) {
      gsl_matrix_set(step->a_mat, i, j, -h * gsl_matrix_get(dfdy, i, j));
    }
    gsl_matrix_set(step->a_mat, i, i, gsl_matrix_get(step->a_mat, i, i) + 1.0);
  }

  /* LU decomposition for the linear system. */
  gsl_linalg_LU_decomp(step->a_mat, step->p_vec, &signum);


  /* Initial step. */

  for(i=0; i<dim; i++) {
    ytemp[i] = h * (yp[i] + h * dfdt[i]);
  }

  gsl_linalg_LU_solve(step->a_mat, step->p_vec, &ytemp_vec, step->delta_temp_vec);

  for (i=0; i<dim; i++) {
    const double di = gsl_vector_get(step->delta_temp_vec, i);
    step->delta[i] = di;
    ytemp[i] = y[i] + di;
  }


  /* Intermediate steps. */

  GSL_ODEIV_FN_EVAL(sys, t, ytemp, y_out);

  for (n_inter=1; n_inter<n_step; n_inter++) {

    for (i=0; i<dim; i++) {
      gsl_vector_set(step->rhs_temp_vec, i, h * y_out[i] - step->delta[i]);
    }

    gsl_linalg_LU_solve(step->a_mat, step->p_vec, step->rhs_temp_vec, step->delta_temp_vec);

    for(i=0; i<dim; i++) {
      step->delta[i] += 2.0 * gsl_vector_get(step->delta_temp_vec, i);
      ytemp[i] += step->delta[i];
    }

    t += h;

    GSL_ODEIV_FN_EVAL(sys, t, ytemp, y_out);
  }


  /* Final step. */

  for(i=0; i<dim; i++) {
    gsl_vector_set(step->rhs_temp_vec, i, h * y_out[i] - step->delta[i]);
  }

  gsl_linalg_LU_solve(step->a_mat, step->p_vec, step->rhs_temp_vec, step->delta_temp_vec);

  for(i=0; i<dim; i++) {
    y_out[i] = ytemp[i] + gsl_vector_get(step->delta_temp_vec, i);
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
    my->k_choice = bsimp_deuf_kchoice(my->eps, my->parent.dimension);
    my->parent.order = 2 * my->k_choice;
    my->h_next   = -GSL_SQRT_DBL_MAX;
  }

  memcpy(my->y_extrap_save, y, my->parent.dimension * sizeof(double));

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
  for(k=0; k <= my->k_choice; k++) {

    const double x_k = (h/bd_sequence[k]) * (h/bd_sequence[k]);

    bsimp_step_local(self,
  		     my->y_extrap_save,
		     my->yp,
  		     my->dfdt,
  		     my->dfdy,
  		     my->parent.dimension,
  		     t_local,
  		     h,
  		     bd_sequence[k],
  		     my->y_extrap_sequence,
  		     sys);

    my->x[k] = x_k;
    poly_extrap(my->d, my->x, k, x_k, my->y_extrap_sequence, y, yerr, my->extrap_work, my->parent.dimension);
  }

  return GSL_SUCCESS;
}
