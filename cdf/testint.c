/* QAGIL: Evaluate an integral over an infinite range using the
   transformation,
   
   integrate(f(x),-Inf,b) = integrate(f(b-(1-t)/t)/t^2,0,1)

   */

struct il_params { double b ; gsl_function * f ; } ;

static double il_transform (double t, void *params);

int
qagil (gsl_function * f,
                       double b,
                       double epsabs, double epsrel, size_t limit,
                       gsl_integration_workspace * workspace,
                       double *result, double *abserr)
{
  int status;

  gsl_function f_transform;
  struct il_params transform_params  ;

  transform_params.b = b ;
  transform_params.f = f ;

  f_transform.function = &il_transform;
  f_transform.params = &transform_params;

  status = gsl_integration_qag (&f_transform, 0.0, 1.0, 
                                epsabs, epsrel, limit, 3,
                                workspace,
                                result, abserr);

  return status;
}

static double 
il_transform (double t, void *params)
{
  struct il_params *p = (struct il_params *) params;
  double b = p->b;
  gsl_function * f = p->f;
  double x = b - (1 - t) / t;
  double y = GSL_FN_EVAL (f, x);
  return (y / t) / t;
}

/* QAGIU: Evaluate an integral over an infinite range using the
   transformation

   integrate(f(x),a,Inf) = integrate(f(a+(1-t)/t)/t^2,0,1)

   */

struct iu_params { double a ; gsl_function * f ; } ;

static double iu_transform (double t, void *params);

int
qagiu (gsl_function * f,
                       double a,
                       double epsabs, double epsrel, size_t limit,
                       gsl_integration_workspace * workspace,
                       double *result, double *abserr)
{
  int status;

  gsl_function f_transform;
  struct iu_params transform_params  ;

  transform_params.a = a ;
  transform_params.f = f ;

  f_transform.function = &iu_transform;
  f_transform.params = &transform_params;

  status = gsl_integration_qag (&f_transform, 0.0, 1.0, 
                                epsabs, epsrel, limit, 3,
                                workspace,
                                result, abserr);

  return status;
}

static double 
iu_transform (double t, void *params)
{
  struct iu_params *p = (struct iu_params *) params;
  double a = p->a;
  gsl_function * f = p->f;
  double x = a + (1 - t) / t;
  double y = GSL_FN_EVAL (f, x);
  return (y / t) / t;
}
