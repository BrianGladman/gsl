/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef _TRIDIAG_H_
#define _TRIDIAG_H_


int solve_tridiag(const double diag[], const double offdiag[], const double b[],
                  double * x,
                  size_t N
                  );

int solve_cyctridiag(const double diag[], const double offdiag[], const double b[],
                     double * x,
                     size_t N
                     );


#endif  /* !_TRIDIAG_H_ */
