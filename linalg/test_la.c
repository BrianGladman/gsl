/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_test.h>
#include <gsl_math.h>
#include "gsl_linalg.h"


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


int
test_matmult(void)
{
  int s = 0;

  gsl_matrix * A = gsl_matrix_calloc(2, 2);
  gsl_matrix * B = gsl_matrix_calloc(2, 3);
  gsl_matrix * C = gsl_matrix_calloc(2, 3);

  gsl_matrix_set(A, 0, 0, 10.0);
  gsl_matrix_set(A, 0, 1,  5.0);
  gsl_matrix_set(A, 1, 0,  1.0);
  gsl_matrix_set(A, 1, 1, 20.0);

  gsl_matrix_set(B, 0, 0, 10.0);
  gsl_matrix_set(B, 0, 1,  5.0);
  gsl_matrix_set(B, 0, 2,  2.0);
  gsl_matrix_set(B, 1, 0,  1.0);
  gsl_matrix_set(B, 1, 1,  3.0);
  gsl_matrix_set(B, 1, 2,  2.0);

  gsl_la_matmult_impl(A, B, C);

  s += ( fabs(gsl_matrix_get(C, 0, 0) - 105.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 0, 1) -  65.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 0, 2) -  30.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 0) -  30.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 1) -  65.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 2) -  42.0) > GSL_DBL_EPSILON );

  return s;
}


int
test_matmult_mod(void)
{
  int s = 0;

  gsl_matrix * A = gsl_matrix_calloc(3, 3);
  gsl_matrix * B = gsl_matrix_calloc(3, 3);
  gsl_matrix * C = gsl_matrix_calloc(3, 3);

  gsl_matrix_set(A, 0, 0, 10.0);
  gsl_matrix_set(A, 0, 1,  5.0);
  gsl_matrix_set(A, 0, 2,  1.0);
  gsl_matrix_set(A, 1, 0,  1.0);
  gsl_matrix_set(A, 1, 1, 20.0);
  gsl_matrix_set(A, 1, 2,  5.0);
  gsl_matrix_set(A, 2, 0,  1.0);
  gsl_matrix_set(A, 2, 1,  3.0);
  gsl_matrix_set(A, 2, 2,  7.0);

  gsl_matrix_set(B, 0, 0, 10.0);
  gsl_matrix_set(B, 0, 1,  5.0);
  gsl_matrix_set(B, 0, 2,  2.0);
  gsl_matrix_set(B, 1, 0,  1.0);
  gsl_matrix_set(B, 1, 1,  3.0);
  gsl_matrix_set(B, 1, 2,  2.0);
  gsl_matrix_set(B, 2, 0,  1.0);
  gsl_matrix_set(B, 2, 1,  3.0);
  gsl_matrix_set(B, 2, 2,  2.0);

  gsl_la_matmult_mod_impl(A, GSL_LA_MOD_NONE, B, GSL_LA_MOD_NONE, C);
  s += ( fabs(gsl_matrix_get(C, 0, 0) - 106.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 0, 1) -  68.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 0, 2) -  32.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 0) -  35.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 1) -  80.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 2) -  52.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 0) -  20.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 1) -  35.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 2) -  22.0) > GSL_DBL_EPSILON );

  gsl_la_matmult_mod_impl(A, GSL_LA_MOD_TRANSPOSE, B, GSL_LA_MOD_NONE, C);
  s += ( fabs(gsl_matrix_get(C, 0, 0) - 102.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 0, 1) -  56.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 0, 2) -  24.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 0) -  73.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 1) -  94.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 2) -  56.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 0) -  22.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 1) -  41.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 2) -  26.0) > GSL_DBL_EPSILON );

  gsl_la_matmult_mod_impl(A, GSL_LA_MOD_NONE, B, GSL_LA_MOD_TRANSPOSE, C);
  s += ( fabs(gsl_matrix_get(C, 0, 0) - 127.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 0, 1) -  27.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 0, 2) -  27.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 0) - 120.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 1) -  71.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 2) -  71.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 0) -  39.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 1) -  24.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 2) -  24.0) > GSL_DBL_EPSILON );

  gsl_la_matmult_mod_impl(A, GSL_LA_MOD_TRANSPOSE, B, GSL_LA_MOD_TRANSPOSE, C);
  s += ( fabs(gsl_matrix_get(C, 0, 0) - 107.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 0, 1) -  15.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 0, 2) -  15.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 0) - 156.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 1) -  71.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 1, 2) -  71.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 0) -  49.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 1) -  30.0) > GSL_DBL_EPSILON );
  s += ( fabs(gsl_matrix_get(C, 2, 2) -  30.0) > GSL_DBL_EPSILON );

  return s;
}


static int
test_LU_solve_dim(size_t dim, const double * actual, double eps)
{
  int s = 0;
  int signum;
  size_t i;

  gsl_vector_int * perm = gsl_vector_int_alloc(dim);
  gsl_vector * rhs = gsl_vector_alloc(dim);
  gsl_matrix * hm  = create_hilbert_matrix(dim);
  gsl_vector * solution = gsl_vector_alloc(dim);
  for(i=0; i<dim; i++) gsl_vector_set(rhs, i, i+1.0);
  s += gsl_la_decomp_LU_impl(hm, perm, &signum);
  s += gsl_la_solve_LU_impl(hm, perm, rhs, solution);
  for(i=0; i<dim; i++) {
    int foo = ( fabs(gsl_vector_get(solution, i) - actual[i])/fabs(actual[i]) > eps );
    if(foo) {
      printf("%3d[%d]: %22.18g   %22.18g\n", dim, i, gsl_vector_get(solution, i), actual[i]);
    }
    s += foo;
  }
  gsl_vector_free(solution);
  gsl_matrix_free(hm);
  gsl_vector_free(rhs);
  gsl_vector_int_free(perm);

  return s;
}


int test_LU_solve(void)
{
  int f;
  int s = 0;
  double actual[16];

  actual[0] =  -8.0;
  actual[1] =  18.0;
  f = test_LU_solve_dim(2, actual, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_LU dim=2");
  s += f;

  actual[0] =   27.0;
  actual[1] = -192.0;
  actual[2] =  210.0;
  f = test_LU_solve_dim(3, actual, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_LU dim=3");
  s += f;

  actual[0] =   -64.0;
  actual[1] =   900.0;
  actual[2] = -2520.0;
  actual[3] =  1820.0;
  f = test_LU_solve_dim(4, actual, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_LU dim=4");
  s += f;

  actual[0]  = -1728.0;
  actual[1]  =  245388.0;
  actual[2]  = -8528520.0;
  actual[3]  =  127026900.0;
  actual[4]  = -1009008000.0;
  actual[5]  =  4768571808.0;
  actual[6]  = -14202796608.0;
  actual[7]  =  27336497760.0;
  actual[8]  = -33921201600.0;
  actual[9]  =  26189163000.0;
  actual[10] = -11437874448.0;
  actual[11] =  2157916488.0;
  f = test_LU_solve_dim(12, actual, 0.05);
  gsl_test(f, "  solve_LU dim=12");
  s += f;

  return s;
}


static int
test_HH_solve_dim(size_t dim, const double * actual, double eps)
{
  int s = 0;
  size_t i;

  gsl_vector * rhs = gsl_vector_alloc(dim);
  gsl_matrix * hm  = create_hilbert_matrix(dim);
  for(i=0; i<dim; i++) gsl_vector_set(rhs, i, i+1.0);
  s += gsl_la_solve_HH_impl(hm, rhs);
  for(i=0; i<dim; i++) {
    int foo = ( fabs(gsl_vector_get(rhs, i) - actual[i])/fabs(actual[i]) > eps );
    if(foo) {
      printf("%3d[%d]: %22.18g   %22.18g\n", dim, i, gsl_vector_get(rhs, i), actual[i]);
    }
    s += foo;
  }
  gsl_matrix_free(hm);
  gsl_vector_free(rhs);

  return s;
}


int test_HH_solve(void)
{
  int f;
  int s = 0;
  double actual[16];

  actual[0] =  -8.0;
  actual[1] =  18.0;
  f = test_HH_solve_dim(2, actual, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_HH dim=2");
  s += f;

  actual[0] =   27.0;
  actual[1] = -192.0;
  actual[2] =  210.0;
  f = test_HH_solve_dim(3, actual, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_HH dim=3");
  s += f;

  actual[0] =   -64.0;
  actual[1] =   900.0;
  actual[2] = -2520.0;
  actual[3] =  1820.0;
  f = test_HH_solve_dim(4, actual, 2.0 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_HH dim=4");
  s += f;

  actual[0]  = -1728.0;
  actual[1]  =  245388.0;
  actual[2]  = -8528520.0;
  actual[3]  =  127026900.0;
  actual[4]  = -1009008000.0;
  actual[5]  =  4768571808.0;
  actual[6]  = -14202796608.0;
  actual[7]  =  27336497760.0;
  actual[8]  = -33921201600.0;
  actual[9]  =  26189163000.0;
  actual[10] = -11437874448.0;
  actual[11] =  2157916488.0;
  f = test_HH_solve_dim(12, actual, 0.01);
  gsl_test(f, "  solve_HH dim=12");
  s += f;

  return s;
}


static int
test_TD_solve_dim(size_t dim, double d, double od, const double * actual, double eps)
{
  int s = 0;
  size_t i;

  gsl_vector * offdiag = gsl_vector_alloc(dim-1);
  gsl_vector * diag = gsl_vector_alloc(dim);
  gsl_vector * rhs = gsl_vector_alloc(dim);
  gsl_vector * sol = gsl_vector_alloc(dim);

  for(i=0; i<dim; i++) {
    gsl_vector_set(diag, i, d);
    gsl_vector_set(rhs,  i, i + 1.0);
  }
  for(i=0; i<dim-1; i++) {
    gsl_vector_set(offdiag, i, od);
  }

  s += gsl_la_solve_symm_tridiag_impl(diag, offdiag, rhs, sol);

  for(i=0; i<dim; i++) {
    double si = gsl_vector_get(sol, i);
    double ai = actual[i];
    int foo = ( fabs(si - ai) / (fabs(ai) + fabs(si)) > eps );
    if(foo) {
      printf("%3d[%d]: %22.18g   %22.18g\n", dim, i, gsl_vector_get(sol, i), actual[i]);
    }
    s += foo;
  }

  gsl_vector_free(sol);
  gsl_vector_free(rhs);
  gsl_vector_free(diag);
  gsl_vector_free(offdiag);

  return s;
}


int test_TD_solve(void)
{
  int f;
  int s = 0;
  double actual[16];

  actual[0] =  0.0;
  actual[1] =  2.0;
  f = test_TD_solve_dim(2, 1.0, 0.5, actual, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_TD dim=2 A");
  s += f;

  actual[0] =  3.0/8.0;
  actual[1] =  15.0/8.0;
  f = test_TD_solve_dim(2, 1.0, 1.0/3.0, actual, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_TD dim=2 B");
  s += f;

  actual[0] =  5.0/8.0;
  actual[1] =  9.0/8.0;
  actual[2] =  2.0;
  actual[3] =  15.0/8.0;
  actual[4] =  35.0/8.0;
  f = test_TD_solve_dim(5, 1.0, 1.0/3.0, actual, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_TD dim=5");
  s += f;

  return s;
}





int main()
{
  gsl_test(test_matmult(),        "Matrix Multiply");
  gsl_test(test_matmult_mod(),    "Matrix Multiply with Modification");
  gsl_test(test_LU_solve(),       "LU Decomposition and Solve");
  gsl_test(test_HH_solve(),       "Householder solve");
  gsl_test(test_TD_solve(),       "Tridiagonal solve");

  return gsl_test_summary();
}
