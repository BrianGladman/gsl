int within_fuzz (double x, double y); /* approximate comparison function */

int 
within_fuzz (const double observed, const double expected)
{
  const double rel = fabs ((observed - expected)/expected) ;
  return (rel < 1e-10);
}
