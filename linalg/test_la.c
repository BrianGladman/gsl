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

gsl_matrix *
create_vandermonde_matrix(int size)
{
  int i, j;
  gsl_matrix * m = gsl_matrix_alloc(size, size);
  for(i=0; i<size; i++) {
    for(j=0; j<size; j++) {
      gsl_matrix_set(m, i, j, pow(i + 1.0, size - j - 1));
    }
  }
  return m;
}

gsl_matrix * hilb2;
gsl_matrix * hilb3;
gsl_matrix * hilb4;
gsl_matrix * hilb12;

double hilb2_solution[] = {-8.0, 18.0} ;
double hilb3_solution[] = {27.0, -192.0, 210.0};
double hilb4_solution[] = {-64.0, 900.0, -2520.0, 1820.0};
double hilb12_solution[] = {-1728.0, 245388.0, -8528520.0, 
                            127026900.0, -1009008000.0, 4768571808.0, 
                            -14202796608.0, 27336497760.0, -33921201600.0,
                            26189163000.0, -11437874448.0, 2157916488.0 };


gsl_matrix * vander2;
gsl_matrix * vander3;
gsl_matrix * vander4;
gsl_matrix * vander12;

double vander2_solution[] = {1.0, 0.0}; 
double vander3_solution[] = {0.0, 1.0, 0.0}; 
double vander4_solution[] = {0.0, 0.0, 1.0, 0.0}; 
double vander12_solution[] = {0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0, 
                            0.0, 0.0, 1.0, 0.0}; 

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

  gsl_matrix_free(A);
  gsl_matrix_free(B);
  gsl_matrix_free(C);

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

  gsl_matrix_free(A);
  gsl_matrix_free(B);
  gsl_matrix_free(C);

  return s;
}


static int
test_LU_solve_dim(const gsl_matrix * m, const double * actual, double eps)
{
  int s = 0;
  int signum;
  size_t i, dim = m->size1;

  gsl_vector_int * perm = gsl_vector_int_alloc(dim);
  gsl_vector * rhs = gsl_vector_alloc(dim);
  gsl_matrix * lu  = gsl_matrix_alloc(dim,dim);
  gsl_vector * solution = gsl_vector_alloc(dim);
  gsl_matrix_copy(lu,m);
  for(i=0; i<dim; i++) gsl_vector_set(rhs, i, i+1.0);
  s += gsl_la_decomp_LU_impl(lu, perm, &signum);
  s += gsl_la_solve_LU_impl(lu, perm, rhs, solution);
  for(i=0; i<dim; i++) {
    int foo = ( fabs(gsl_vector_get(solution, i) - actual[i])/fabs(actual[i]) > eps );
    if(foo) {
      printf("%3d[%d]: %22.18g   %22.18g\n", dim, i, gsl_vector_get(solution, i), actual[i]);
    }
    s += foo;
  }
  gsl_vector_free(solution);
  gsl_matrix_free(lu);
  gsl_vector_free(rhs);
  gsl_vector_int_free(perm);

  return s;
}


int test_LU_solve(void)
{
  int f;
  int s = 0;

  f = test_LU_solve_dim(hilb2, hilb2_solution, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_LU hilbert(2)");
  s += f;

  f = test_LU_solve_dim(hilb3, hilb3_solution, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_LU hilbert(3)");
  s += f;

  f = test_LU_solve_dim(hilb4, hilb4_solution, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_LU hilbert(4)");
  s += f;

  f = test_LU_solve_dim(hilb12, hilb12_solution, 0.05);
  gsl_test(f, "  solve_LU hilbert(12)");
  s += f;

  f = test_LU_solve_dim(vander2, vander2_solution, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_LU vander(2)");
  s += f;

  f = test_LU_solve_dim(vander3, vander3_solution, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_LU vander(3)");
  s += f;

  f = test_LU_solve_dim(vander4, vander4_solution, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_LU vander(4)");
  s += f;

  f = test_LU_solve_dim(vander12, vander12_solution, 0.05);
  gsl_test(f, "  solve_LU vander(12)");
  s += f;

  return s;
}



static int
test_QR_solve_dim(const gsl_matrix * m, const double * actual, double eps)
{
  int s = 0;
  int signum;
  size_t i, dim = m->size1;

  gsl_vector_int * perm = gsl_vector_int_alloc(dim);
  gsl_vector * rhs = gsl_vector_alloc(dim);
  gsl_matrix * qr  = gsl_matrix_alloc(dim,dim);
  gsl_vector * d = gsl_vector_alloc(dim);
  gsl_vector * solution = gsl_vector_alloc(dim);
  gsl_matrix_copy(qr,m);
  for(i=0; i<dim; i++) gsl_vector_set(rhs, i, i+1.0);
  s += gsl_la_decomp_QR_impl(qr, d, perm, &signum);
  s += gsl_la_solve_QR_impl(qr, d, perm, rhs, solution);
  for(i=0; i<dim; i++) {
    int foo = ( fabs(gsl_vector_get(solution, i) - actual[i])/fabs(actual[i]) > eps );
    if(foo) {
      printf("%3d[%d]: %22.18g   %22.18g\n", dim, i, gsl_vector_get(solution, i), actual[i]);
    }
    s += foo;
  }
  gsl_vector_free(solution);
  gsl_vector_free(d);
  gsl_matrix_free(qr);
  gsl_vector_free(rhs);
  gsl_vector_int_free(perm);

  return s;
}

int test_QR_solve(void)
{
  int f;
  int s = 0;

  f = test_QR_solve_dim(hilb2, hilb2_solution, 2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_QR hilbert(2)");
  s += f;

  f = test_QR_solve_dim(hilb3, hilb3_solution, 2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_QR hilbert(3)");
  s += f;

  f = test_QR_solve_dim(hilb4, hilb4_solution, 2 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_QR hilbert(4)");
  s += f;

  f = test_QR_solve_dim(hilb12, hilb12_solution, 0.05);
  gsl_test(f, "  solve_QR hilbert(12)");
  s += f;

  f = test_QR_solve_dim(vander2, vander2_solution, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_QR vander(2)");
  s += f;

  f = test_QR_solve_dim(vander3, vander3_solution, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_QR vander(3)");
  s += f;

  f = test_QR_solve_dim(vander4, vander4_solution, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_QR vander(4)");
  s += f;

  f = test_QR_solve_dim(vander12, vander12_solution, 0.05);
  gsl_test(f, "  solve_QR vander(12)");
  s += f;

  return s;
}


static int
test_QR_update_dim(const gsl_matrix * m, double eps)
{
  int s = 0;
  int signum;
  size_t i,j, dim = m->size1;

  gsl_vector_int * perm = gsl_vector_int_alloc(dim);
  gsl_vector * rhs = gsl_vector_alloc(dim);
  gsl_matrix * qr1  = gsl_matrix_alloc(dim,dim);
  gsl_matrix * qr2  = gsl_matrix_alloc(dim,dim);
  gsl_matrix * q  = gsl_matrix_alloc(dim,dim);
  gsl_matrix * r  = gsl_matrix_alloc(dim,dim);
  gsl_vector * d = gsl_vector_alloc(dim);
  gsl_vector * solution1 = gsl_vector_alloc(dim);
  gsl_vector * solution2 = gsl_vector_alloc(dim);
  gsl_vector * u = gsl_vector_alloc(dim);
  gsl_vector * v = gsl_vector_alloc(dim);
  gsl_vector * w = gsl_vector_alloc(dim);
  gsl_matrix_copy(qr1,m);
  gsl_matrix_copy(qr2,m);
  for(i=0; i<dim; i++) gsl_vector_set(rhs, i, i+1.0);
  for(i=0; i<dim; i++) gsl_vector_set(u, i, sin(i+1.0));
  for(i=0; i<dim; i++) gsl_vector_set(v, i, sin(i*i+3.0));

  for(i=0; i<dim; i++) 
    {
      double ui = gsl_vector_get(u, i);
      for(j=0; j<dim; j++) 
        {
          double vj = gsl_vector_get(v, j);
          double qij = gsl_matrix_get(qr1, i, j);
          gsl_matrix_set(qr1, i, j, qij + ui * vj);
        }
    }

  s += gsl_la_decomp_QR_impl(qr1, d, perm, &signum);
  s += gsl_la_solve_QR_impl(qr1, d, perm, rhs, solution1);

  s += gsl_la_decomp_QR_impl(qr2, d, perm, &signum);
  s += gsl_la_unpack_QR_impl(qr2, d, q, r);
  s += gsl_la_update_QR_impl(q, r, perm, u, v, w);
  s += gsl_la_qrsolve_QR_impl(q, r, perm, rhs, solution2);

  for(i=0; i<dim; i++) {
    double s1 = gsl_vector_get(solution1, i);
    double s2 = gsl_vector_get(solution2, i);

    int foo = ( fabs(s1 - s2)/fabs(s2) > eps );
    if(foo) {
      printf("%3d[%d]: %22.18g   %22.18g\n", dim, i, s1, s2);
    }
    s += foo;
  }
  gsl_vector_free(solution1);
  gsl_vector_free(solution2);
  gsl_vector_free(d);
  gsl_vector_free(u);
  gsl_vector_free(v);
  gsl_vector_free(w);
  gsl_matrix_free(qr1);
  gsl_matrix_free(qr2);
  gsl_matrix_free(q);
  gsl_matrix_free(r);
  gsl_vector_free(rhs);
  gsl_vector_int_free(perm);

  return s;
}

int test_QR_update(void)
{
  int f;
  int s = 0;

  f = test_QR_update_dim(hilb2,  2 * 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  update_QR hilbert(2)");
  s += f;

  f = test_QR_update_dim(hilb3,  2 * 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  update_QR hilbert(3)");
  s += f;

  f = test_QR_update_dim(hilb4, 2 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  update_QR hilbert(4)");
  s += f;

  f = test_QR_update_dim(hilb12, 0.05);
  gsl_test(f, "  update_QR hilbert(12)");
  s += f;

  f = test_QR_update_dim(vander2, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  update_QR vander(2)");
  s += f;

  f = test_QR_update_dim(vander3, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  update_QR vander(3)");
  s += f;

  f = test_QR_update_dim(vander4, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  update_QR vander(4)");
  s += f;

  f = test_QR_update_dim(vander12, 0.05);
  gsl_test(f, "  update_QR vander(12)");
  s += f;

  return s;
}


static int
test_HH_solve_dim(const gsl_matrix * m, const double * actual, double eps)
{
  int s = 0;
  size_t i, dim = m->size1;

  gsl_vector_int * perm = gsl_vector_int_alloc(dim);
  gsl_matrix * hh  = gsl_matrix_alloc(dim,dim);
  gsl_vector * d = gsl_vector_alloc(dim);
  gsl_vector * solution = gsl_vector_alloc(dim);
  gsl_matrix_copy(hh,m);
  for(i=0; i<dim; i++) gsl_vector_set(solution, i, i+1.0);
  s += gsl_la_solve_HH_impl(hh, solution);
  for(i=0; i<dim; i++) {
    int foo = ( fabs(gsl_vector_get(solution, i) - actual[i])/fabs(actual[i]) > eps );
    if( foo) {
      printf("%3d[%d]: %22.18g   %22.18g\n", dim, i, gsl_vector_get(solution, i), actual[i]);
    }
    s += foo;
  }
  gsl_vector_free(solution);
  gsl_vector_free(d);
  gsl_matrix_free(hh);
  gsl_vector_int_free(perm);

  return s;
}

int test_HH_solve(void)
{
  int f;
  int s = 0;

  f = test_HH_solve_dim(hilb2, hilb2_solution, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_HH hilbert(2)");
  s += f;

  f = test_HH_solve_dim(hilb3, hilb3_solution, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_HH hilbert(3)");
  s += f;

  f = test_HH_solve_dim(hilb4, hilb4_solution, 2.0 * 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_HH hilbert(4)");
  s += f;

  f = test_HH_solve_dim(hilb12, hilb12_solution, 0.01);
  gsl_test(f, "  solve_HH hilbert(12)");
  s += f;

  f = test_HH_solve_dim(vander2, vander2_solution, 8.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_HH vander(2)");
  s += f;

  f = test_HH_solve_dim(vander3, vander3_solution, 64.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_HH vander(3)");
  s += f;

  f = test_HH_solve_dim(vander4, vander4_solution, 1024.0 * GSL_DBL_EPSILON);
  gsl_test(f, "  solve_HH vander(4)");
  s += f;

  f = test_HH_solve_dim(vander12, vander12_solution, 0.05);
  gsl_test(f, "  solve_HH vander(12)");
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
  hilb2 = create_hilbert_matrix(2);
  hilb3 = create_hilbert_matrix(3);
  hilb4 = create_hilbert_matrix(4);
  hilb12 = create_hilbert_matrix(12);

  vander2 = create_vandermonde_matrix(2);
  vander3 = create_vandermonde_matrix(3);
  vander4 = create_vandermonde_matrix(4);
  vander12 = create_vandermonde_matrix(12);

  gsl_test(test_matmult(),        "Matrix Multiply");
  gsl_test(test_matmult_mod(),    "Matrix Multiply with Modification");
  gsl_test(test_LU_solve(),       "LU Decomposition and Solve");
  gsl_test(test_QR_solve(),       "QR Decomposition and Solve");
  gsl_test(test_QR_update(),      "QR Rank-1 Update");
  gsl_test(test_HH_solve(),       "Householder solve");
  gsl_test(test_TD_solve(),       "Tridiagonal solve");

  gsl_matrix_free(hilb2);
  gsl_matrix_free(hilb3);
  gsl_matrix_free(hilb4);
  gsl_matrix_free(hilb12);

  gsl_matrix_free(vander2);
  gsl_matrix_free(vander3);
  gsl_matrix_free(vander4);
  gsl_matrix_free(vander12);

  return gsl_test_summary();
}
