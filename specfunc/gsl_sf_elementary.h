/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Miscellaneous elementary functions and operations.
 */
#ifndef __GSL_SF_ELEMENTARY_H__
#define __GSL_SF_ELEMENTARY_H__

#include <gsl/gsl_sf_result.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


/* Multiplication.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_multiply_impl(const double x, const double y, gsl_sf_result * result);
int gsl_sf_multiply_e(const double x, const double y, gsl_sf_result * result);


/* Multiplication of quantities with associated errors.
 */
int gsl_sf_multiply_err_impl(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);
int gsl_sf_multiply_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_ELEMENTARY_H__ */
