/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef TEST_SF_H
#define TEST_SF_H

#include <gsl_machine.h>
#include "gsl_sf_result.h"

int test_sf_analyze(const char * message, gsl_sf_result r,
                    double expect_val, double tol,
                    int return_val, int expect_return);

#define TEST_TOL0  (2.0*GSL_DBL_EPSILON)
#define TEST_TOL1  (16.0*GSL_DBL_EPSILON)
#define TEST_TOL2  (128.0*GSL_DBL_EPSILON)
#define TEST_TOL3  (1024.0*GSL_DBL_EPSILON)


#define TEST_SF(stat, func, args, val, tol, expect_return)                   \
do {                                                                         \
  int _tsf_status = func args;                                               \
  stat += test_sf_analyze(#func, r, val, tol, _tsf_status, expect_return);   \
} while(0)



int test_airy(void);


#endif /* !TEST_SF_H */
