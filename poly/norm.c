static double norm_companion (const double *m, size_t nc);

static double
norm_companion (const double *m, size_t nc)
{
  size_t i;
  double A, s2 = 0;

  for (i = 1; i < nc - 1; i++)
    {
      const double t = MAT (m, i + 1, i, nc);
      s2 += t * t;
    }

  for (i = 0; i < nc; i++)
    {
      const double t = MAT (m, i, nc, nc);
      s2 += t * t;
    }

  A = sqrt (s2);

  return A;
}

