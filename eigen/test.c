/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include "gsl_eigen.h"

gsl_matrix *
create_hilbert_matrix(int size)
{
  int i, j;
  gsl_matrix * m = gsl_matrix_alloc(size, size);
  for(i=0; i<size; i++) {
    for(j=0; j<size; j++) {
      gsl_matrix_set(m, i, j, 1.0/(i+j+1.0));
    }
  }
  return m;
}

int test_eigen_jacobi(void)
{
  int s = 0;

  int nrot;
  gsl_matrix * evec = gsl_matrix_alloc(10, 10);
  gsl_vector * eval = gsl_vector_alloc(10);

  /* eigenvalues of 10x10 Hilbert matrix */
  const double eval_h10[10] = { 1.7519196702651775224,
                                0.3429295484835090962,
                                0.03574181627163923589,
                                0.0025308907686700381437,
                                0.00012874961427637707981,
                                4.7296892931823475060e-06,
                                1.2289677387511750496e-07,
                                2.1474388173504786077e-09,
                                2.2667467477629255254e-11,
                                1.0931538193796657185e-13
                              };

  /* 10x10 Hilbert matrix */
  gsl_matrix * hm = create_hilbert_matrix(10);
  gsl_eigen_jacobi_impl(hm, eval, evec, 1000, &nrot);
  gsl_eigen_sort_impl(eval, evec, GSL_EIGEN_SORT_VALUE);

  s += ( fabs(eval_h10[0] - eval->data[9]) > 1.0e-15 );
  s += ( fabs(eval_h10[1] - eval->data[8]) > 1.0e-14 );
  s += ( fabs(eval_h10[2] - eval->data[7]) > 1.0e-13 );
  /* */
  s += ( fabs(eval_h10[8] - eval->data[1]) > 1.0e-04 );
  s += ( fabs(eval_h10[9] - eval->data[0]) > 1.0e-03 );

/* FIXME: must check eigenvectors as well */

  gsl_matrix_free(hm);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);

  return s;
}

int test_invert_jacobi(void)
{
  int s = 0;
  int i, j;
  gsl_matrix * hminv = gsl_matrix_alloc(10, 10);
  gsl_matrix * id    = gsl_matrix_alloc(10, 10);

  /* 10x10 Hilbert matrix */
  gsl_matrix * hm = create_hilbert_matrix(10);
  gsl_eigen_invert_jacobi_impl(hm, hminv, 1000);

  gsl_la_matmult_impl(hm, hminv, id);

  for(i=0; i<10; i++) {
    for(j=0; j<10; j++) {
      double delta_ij = ( i == j ? 1.0 : 0.0 );
      double id_ij    = gsl_matrix_get(id, i, j);
      int rs = ( fabs(id_ij - delta_ij) > 5.0e-3 );
      s += rs;
    }
  }

  gsl_matrix_free(hm);
  gsl_matrix_free(hminv);
  gsl_matrix_free(id);

  return s;
}


int main()
{
  gsl_test(test_eigen_jacobi(),   "Eigensystem:  Jacobi Method");
  gsl_test(test_invert_jacobi(),  "Inversion:    Jacobi Method");

  return gsl_test_summary();
}
