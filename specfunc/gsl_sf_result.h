/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_RESULT_H_
#define GSL_SF_RESULT_H_

struct gsl_sf_result_struct {
  double val;
  double err;
  int    _reserved;
};
typedef struct gsl_sf_result_struct gsl_sf_result;


#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


#endif /* GSL_SF_RESULT_H_ */
