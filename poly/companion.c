static void set_companion_matrix (const double *a, size_t n, double *m);

static void
set_companion_matrix (const double *a, size_t nc, double *m)
{
  size_t i, j;

  for (i = 0; i < nc; i++)
    for (j = 0; j < nc; j++)
      MAT (m, i, j, nc) = 0.0;

  for (i = 1; i < nc; i++)
    MAT (m, i, i - 1, nc) = 1.0;

  for (i = 0; i < nc; i++)
    MAT (m, i, nc - 1, nc) = -a[i] / a[nc];
}
