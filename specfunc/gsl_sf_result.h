/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_RESULT_H_
#define GSL_SF_RESULT_H_

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_impl(const gsl_sf_result_e10 * re, gsl_sf_result * r);
int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


#endif /* GSL_SF_RESULT_H_ */
