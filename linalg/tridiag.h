/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Low level tridiagonal solvers.
 * Used internally in other areas of GSL.
 */
#ifndef _GSL_TRIDIAG_H_
#define _GSL_TRIDIAG_H_

#include <stdlib.h>


int solve_tridiag(
  const double diag[], size_t d_stride,
  const double offdiag[], size_t o_stride,
  const double b[], size_t b_stride,
  double x[], size_t x_stride,
  size_t N
  );

int solve_cyctridiag(
  const double diag[], size_t d_stride,
  const double offdiag[], size_t o_stride,
  const double b[], size_t b_stride,
  double x[], size_t x_stride,
  size_t N
  );


#endif  /* !_GSL_TRIDIAG_H_ */
