int within_fuzz (double x, double y); /* approximate comparison function */

int 
within_fuzz (double observed, double expected)
{
  return fabs ((observed - expected)/expected) < 0.00001;
}
