/* Generate a Givens rotation (cos,sin) which takes v=(x,y) to (|v|,0) 

   From Golub and Van Loan, "Matrix Computations", Section 5.1.8 */

inline static void
create_givens (const double a, const double b, double *c, double *s)
{
  if (b == 0)
    {
      *c = 1;
      *s = 0;
    }
  else if (fabs (b) > fabs (a))
    {
      double t = -a / b;
      double s1 = 1.0 / sqrt (1 + t * t);
      *s = s1;
      *c = s1 * t;
    }
  else
    {
      double t = -b / a;
      double c1 = 1.0 / sqrt (1 + t * t);
      *c = c1;
      *s = c1 * t;
    }
}

inline static void
apply_givens_qr (size_t M, size_t N, gsl_matrix * q, gsl_matrix * r,
		 size_t i, size_t j, double c, double s)
{
  size_t k;

  /* Apply rotation to matrix Q,  Q' = Q G */

  for (k = 0; k < N; k++)
    {
      double qki = gsl_matrix_get (q, k, i);
      double qkj = gsl_matrix_get (q, k, j);
      gsl_matrix_set (q, k, i, qki * c - qkj * s);
      gsl_matrix_set (q, k, j, qki * s + qkj * c);
    }

  /* Apply rotation to matrix R, R' = G^T R (note: upper triangular so
     zero for column < row) */

  for (k = GSL_MIN (i, j); k < N; k++)
    {
      double rik = gsl_matrix_get (r, i, k);
      double rjk = gsl_matrix_get (r, j, k);
      gsl_matrix_set (r, i, k, c * rik - s * rjk);
      gsl_matrix_set (r, j, k, s * rik + c * rjk);
    }
}

inline static void
apply_givens_vec (gsl_vector * v, size_t i, size_t j, double c, double s)
{
  /* Apply rotation to vector v' = G^T v */

  double vi = gsl_vector_get (v, i);
  double vj = gsl_vector_get (v, j);
  gsl_vector_set (v, i, c * vi - s * vj);
  gsl_vector_set (v, j, s * vi + c * vj);
}
