
static inline int
cheb_eval_e(const cheb_series * cs,
            const double x,
            gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }

  result->val = y*d - dd + 0.5 * cs->c[0];
  result->err = GSL_DBL_EPSILON * fabs(result->val) + fabs(cs->c[cs->order]);
  return GSL_SUCCESS;
}

