/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_ODEIV_H
#define GSL_ODEIV_H

#include <stdio.h>
#include <stdlib.h>


/* Description of a system of ODEs.
 *
 * y' = f(t,y) = dydt(t, y)
 *
 * The system is specified by giving the right-hand-side
 * of the equation and possibly a jacobian function.
 *
 * Some methods require the jacobian function, which calculates
 * the matrix dfdy and the vector dfdt. The matrix dfdy conforms
 * to the GSL standard, being a continuous range of floating point
 * values, in row-order.
 *
 * As with GSL function objects, user-supplied parameter
 * data is also present. 
 */
typedef struct  {
  int (* function) (double t, const double y[], double dydt[], void * params);
  int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params);
  size_t dimension;
  void * params;
}
gsl_odeiv_system;

#define GSL_ODEIV_FN_EVAL(S,t,y,f)  (*((S)->function))(t,y,f,(S)->params)
#define GSL_ODEIV_JA_EVAL(S,t,y,dfdy,dfdt)  (*((S)->jacobian))(t,y,dfdy,dfdt,(S)->params)


/* General stepper object.
 *
 * Opaque object for stepping an ODE system from t to t+h.
 * In general the object has some state which facilitates
 * iterating the stepping operation.
 */
typedef struct {
  char * _name;
  int  (*_step)  (void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
  int  (*_reset) (void * _state);
  void (*_free)  (void * _state, void * _work);
  void * _state;
  void * _work;
  int can_use_dydt;
  size_t dimension;
  unsigned int order;
  int stutter;
}
gsl_odeiv_step;


/* General stepper object factory.
 *
 * Specialized factory instances create
 * steppers embodying different methods.
 */
typedef struct {
  const char * name;
  gsl_odeiv_step * (*create) (unsigned int dimension);
}
gsl_odeiv_step_factory;


/* Available stepper factories.
 *
 * rk2    : embedded 2nd(3rd) Runge-Kutta
 * rk4    : 4th order (classical) Runge-Kutta
 * rkck   : embedded 4th(5th) Runge-Kutta, Cash-Karp
 * rk8pd  : embedded 8th(9th) Runge-Kutta, Prince-Dormand
 * rk2imp : implicit 2nd order Runge-Kutta at Gaussian points
 * rk4imp : implicit 4th order Runge-Kutta at Gaussian points
 * gear1  : M=1 implicit Gear method
 * gear2  : M=2 implicit Gear method
 */
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_rk2;
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_rk4;
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_rkck;
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_rk8pd;
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_rk2imp;
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_rk4imp;
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_gear1;
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_gear2;


/* General stepper object methods.
 */
const char * gsl_odeiv_step_name(const gsl_odeiv_step * s);
int  gsl_odeiv_step_impl(gsl_odeiv_step * s, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
int  gsl_odeiv_step_reset(gsl_odeiv_step * s);
void gsl_odeiv_step_free(gsl_odeiv_step * s);


/* General evolution monitor object.
 *
 * Instances of this object are used by
 * the evolution to facilitate a kind of
 * callback mechanism.
 */
typedef struct {
  int (*pre_step)  (void * self, double t, unsigned int dim, const double y[]);
  int (*post_step) (void * self, double t, unsigned int dim, const double y[], const double yerr[]);
  FILE * f;
  void * params;
}
gsl_odeiv_evolve_mon;


/* Specialized monitor constructors.
 */
gsl_odeiv_evolve_mon * gsl_odeiv_evolve_mon_stream_new(FILE *);

void gsl_odeiv_evolve_mon_free(gsl_odeiv_evolve_mon *);


/* General step size control object.
 *
 * The hadjust() method controls the adjustment of
 * step size given the result of a step and the error.
 * Valid hadjust() methods must return one of the codes below.
 *
 * The general data can be used by specializations
 * to store state and control their heuristics.
 */
typedef struct {
  int  (*_hadjust) (void * data, size_t dim, unsigned int ord, const double y[], const double yerr[], const double yp[], double * h);
  void (*_free) (void * data);
  size_t _data_size;
  void * _data;
}
gsl_odeiv_evolve_control;


/* Possible return values for an hadjust() evolution method.
 */
#define GSL_ODEIV_HADJ_INC   1  /* step was increased */
#define GSL_ODEIV_HADJ_NIL   0  /* step unchanged     */
#define GSL_ODEIV_HADJ_DEC (-1) /* step decreased     */

#define GSL_ODEIV_CONTROL_HADJ(c,dim,ord,y,yerr,yp,h) (*((c)->_hadjust))((c)->_data,dim,ord,y,yerr,yp,h)


/* Available control object constructors.
 *
 * The standard control object is a four parameter heuristic
 * defined as follows:
 *    D0 = eps_rel * (a_y |y| + a_dydt h |y'|) + eps_abs
 *    D1 = |yerr|
 *    q  = consistency order of method (q=4 for 4(5) embedded RK)
 *    S  = safety factor (0.9 say)
 *
 *                      /  (D0/D1)^(1/(q+1))  D0 >= D1
 *    h_NEW = S h_OLD * |
 *                      \  (D0/D1)^(1/q)      D0 < D1
 *
 * This encompasses all the standard error scaling methods.
 *
 * The y method is the standard method with a_y=1, a_dydt=0.
 * The yp method is the standard method with a_y=0, a_dydt=1.
 */

gsl_odeiv_evolve_control * gsl_odeiv_evolve_control_standard_new(double eps_rel, double eps_abs, double a_y, double a_dydt);
gsl_odeiv_evolve_control * gsl_odeiv_evolve_control_y_new(double eps_rel, double eps_abs);
gsl_odeiv_evolve_control * gsl_odeiv_evolve_control_yp_new(double eps_rel, double eps_abs);

void gsl_odeiv_evolve_control_free(gsl_odeiv_evolve_control *);



/* General evolution object.
 */
typedef struct {
  double * y0;
  double * yerr;
  double * dydt_in;
  double * dydt_out;
  size_t dimension;
  unsigned int count;
  unsigned int count_stutter;
}
gsl_odeiv_evolve;


/* Evolution object methods.
 */
gsl_odeiv_evolve * gsl_odeiv_evolve_new(void);
int gsl_odeiv_evolve_impl(gsl_odeiv_evolve *, gsl_odeiv_evolve_mon * mon, gsl_odeiv_evolve_control * con, gsl_odeiv_step * step, const gsl_odeiv_system * dydt, double t0, double t1, double hstart, double y[]);
void gsl_odeiv_evolve_free(gsl_odeiv_evolve *);


#endif  /* !GSL_ODEIV_H */
