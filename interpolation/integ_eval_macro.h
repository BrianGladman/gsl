/* stupid littel macro for doing the spline integral evaluation
   which is common to both the cspline and akima methods
 */


#define INTEG_EVAL(ai, bi, ci, di, xi, a, b, res)     \
do {                                                  \
  const double _bma = (b)-(a);                        \
  const double _t0 = (b) + (a);                       \
  const double _t1 = (a)*(a) + (a)*(b) + (b)*(b);                             \
  const double _t2 = (a)*(a)*(a) + (a)*(a)*(b) + (b)*(b)*(a) + (b)*(b)*(b);   \
  const double _bterm = 0.5 * (bi) * (_t0-2.0*(xi));                          \
  const double _cterm = (ci) / 3.0 * (_t1 - 3.0*(xi)*(_t0 - (xi)));           \
  const double _dterm = (di) / 4.0 * (_t2 - 2.0*(xi)*(2.0*_t1 - (xi)*(3.0*_t0 - 2.0*(xi))));  \
  res = _bma * ((ai) + _bterm + _cterm + _dterm);  \
} while (0)
