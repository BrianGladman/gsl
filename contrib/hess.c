int
hess (const unsigned int N, gsl_matrix * A)
{
  /*
     Brings a given square matrix A to a similar upper Hessenberg.

     N is the order of A. For sure, there should be a way to
     determine the dimensions of A, but I am not so used yet with
     GSL.

     Ricardo Biloti <biloti@mat.ufpr.br>
     Department of Mathematics
     Federal University of Parana - Brazil
   */

  gsl_matrix_view AA;
  gsl_vector *v_full;
  gsl_vector_view x, v;

  double tau_i;
  unsigned int j;

  v_full = gsl_vector_alloc (N);
  for (j = 0; j < N - 2; j++)
    {
      x = gsl_matrix_column (A, j);
      gsl_vector_memcpy (v_full, &(x.vector));
      v = gsl_vector_subvector (v_full, j + 1, N - j - 1);

      tau_i = gsl_linalg_householder_transform (&(v.vector));

      AA = gsl_matrix_submatrix (A, j + 1, j, N - j - 1, N - j);
      gsl_linalg_householder_hm (tau_i, &(v.vector), &(AA.matrix));

      AA = gsl_matrix_submatrix (A, 0, j + 1, N, N - j - 1);
      gsl_linalg_householder_mh (tau_i, &(v.vector), &(AA.matrix));

    }
  gsl_vector_free (v_full);

  return EXIT_SUCCESS;
}
