/* Author:  G. Jungman */
/*
 * Automatic generation of error handling versions
 * of functions. Given a "func_name", these will
 * create the implementations of
 *      func_name_e()
 *      func_name()
 * in terms of func_name_impl().
 *
 * Example: MAKE_FUNC_ERRHAND(gsl_foo, (double x, double * result), (x, result))
 *          should expand to
 *
 *          int gsl_foo_e(double x, double * result)
 *          {
 *            int status = gsl_foo_impl(a, result);
 *            if(status != GSL_SUCCESS) {
 *              GSL_ERROR("gsl_foo", status);
 *            }
 *            return status;
 *          }
 *
 *          and MAKE_FUNC_NATURAL(gsl_foo, (double x), (x, &y))
 *          should expand to
 *
 *          double gsl_foo(double x)
 *          {
 *            double y;
 *            int status = gsl_foo_impl(x, &y);
 *            if(status != GSL_SUCCESS) {
 *              GSL_WARNING("gsl_foo", status);
 *            }
 *            return y;
 *          }
 *
 */
#ifndef _TEMPLATES_ERRFUNCS_H_
#define _TEMPLATES_ERRFUNCS_H_


#define NAME_IMPL(f)  f ## _ ## impl
#define NAME_E(f)     f ## _ ## e


#define MAKE_FUNC_ERRHAND(func, args, impl_args)     \
                                                     \
int NAME_E(func) args                                \
{                                                    \
  int status = NAME_IMPL(func) impl_args;            \
  if(status != GSL_SUCCESS) {                        \
    GSL_ERROR(#func "_e", status);                   \
  }                                                  \
  return status;                                     \
}                                                    \


#define MAKE_FUNC_NATURAL(func, args, impl_args)     \
                                                     \
double func args                                     \
{                                                    \
  double y;                                          \
  int status = NAME_IMPL(func) impl_args;            \
  if(status != GSL_SUCCESS) {                        \
    GSL_WARNING(#func, status);                      \
  }                                                  \
  return y;                                          \
}                                                    \


#endif /* !_TEMPLATES_ERRFUNCS_H_ */
