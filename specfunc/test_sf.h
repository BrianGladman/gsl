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

#define TEST_SF_INCONS  1
#define TEST_SF_ERRNEG  2
#define TEST_SF_TOLBAD  4
#define TEST_SF_RETBAD  8


#define TEST_SF(stat, func, args, val, tol, expect_return)                        \
do {                                                                              \
  char message_buff[4096];                                                  	  \
  int _tsf_local_s = 0;                                                     	  \
  int _tsf_status = func args;                                              	  \
  _tsf_local_s |= test_sf_check_result(message_buff, r, val, tol);	          \
  _tsf_local_s |= test_sf_check_return(message_buff, _tsf_status, expect_return); \
  gsl_test(_tsf_local_s, #func #args);                                        	  \
  if(_tsf_local_s != 0) {                                                         \
    printf("%s", message_buff);                                                   \
  }                                                                               \
  stat += _tsf_local_s;                                                           \
} while(0)


#define TEST_SF_2(stat, func, args, val1, tol1, val2, tol2, expect_return)  	  \
do {                                                                        	  \
  char message_buff[4096];                                                  	  \
  int _tsf_local_s = 0;                                                     	  \
  int _tsf_status = func args;                                              	  \
  _tsf_local_s |= test_sf_check_result(message_buff, r1, val1, tol1);	          \
  _tsf_local_s |= test_sf_check_result(message_buff, r2, val2, tol2);	          \
  _tsf_local_s |= test_sf_check_return(message_buff, _tsf_status, expect_return); \
  gsl_test(_tsf_local_s, #func #args);                                        	  \
  if(_tsf_local_s != 0) {                                                         \
    printf("%s", message_buff);                                                   \
  }                                                                               \
  stat += _tsf_local_s;                                                           \
} while(0)


#define TEST_SF_SGN(stat, func, args, val_in, tol, expect_sgn, expect_return)     \
do {                                                                              \
  char message_buff[4096];                                                  	  \
  int _tsf_local_s = 0;                                                     	  \
  int _tsf_status = func args;                                              	  \
  gsl_sf_result _tsf_local_r;                                             	  \
  _tsf_local_r.val = sgn;                                                         \
  _tsf_local_r.err = 0.0;                                             	          \
  _tsf_local_s |= test_sf_check_result(message_buff, r, val_in, tol);             \
  _tsf_local_s |= test_sf_check_result(message_buff, _tsf_local_r, expect_sgn, 0.0); \
  _tsf_local_s |= test_sf_check_return(message_buff, _tsf_status, expect_return); \
  gsl_test(_tsf_local_s, #func #args);                                        	  \
  if(_tsf_local_s != 0) {                                                         \
    printf("%s", message_buff);                                                   \
  }                                                                               \
  stat += _tsf_local_s;                                                           \
} while(0)


int test_airy(void);


#endif /* !TEST_SF_H */
