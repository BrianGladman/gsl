/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_ODEIV_H
#define GSL_ODEIV_H


/* function type for right hand side of ODE */
struct gsl_odeiv_function_struct {
  int (* function) (double t, const double y[], double f[], void * params);
  void * params;
};
typedef struct gsl_odeiv_function_struct gsl_odeiv_function;


/* general stepper object */
struct _gsl_odeiv_step_struct {
  int  (*_step)  (void *, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt);
  int  (*_reset) (void *);
  void (*_free)  (void *);
  void * _state;
  unsigned int dimension;
};
typedef struct _gsl_odeiv_step_struct gsl_odeiv_step;

#define GSL_ODEIV_FN_EVAL(F,t,y,f)  (*((F)->function))(t,y,f,(F)->params)


/* stepper object factory */
typedef struct {
  const char * name;
  gsl_odeiv_step * (*create) (unsigned int dimension);
}
gsl_odeiv_step_factory;


/* available stepper factories */
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_rk4;
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_rkck;
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_rk4imp;
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_gear1;
extern const gsl_odeiv_step_factory  gsl_odeiv_step_factory_gear2;


/* stepper methods */

int  gsl_odeiv_step_impl(gsl_odeiv_step * s, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt);
int  gsl_odeiv_step_reset(gsl_odeiv_step * s);
void gsl_odeiv_step_free(gsl_odeiv_step * s);


#endif  /* !GSL_ODEIV_H */
