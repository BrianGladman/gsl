/* Created: [GJ] Sat Apr 27 01:56:04 EDT 1996
 *
 * [GJ] Mon Feb 10 10:16:40 EST 1997
 * Major overhaul of interface.
 */
#ifndef MATRICES_H_
#define MATRICES_H_
#include <stdio.h>
#include <string.h>


/* Symmetric matrix diagonalization
 * doctored from Numerical Recipes,
 * with a sane indexing convention
 * and no required external routines.
 *
 * The columns of v[][] are the returned eigenvectors.
 * d[] is the returned eigenvalues.
 *
 * I have checked this on some simple examples. (GJ)
 *
 * WARNING: jacobi_diag() destroys the super-diagonal elements of the
 *          passed matrix; if this matrix is needed again it will
 *          have to be re-constructed from its lower triangular part
 *
 * The arguments are assumed to be already allocated.
 */
extern void jacobi_diag(double ** m, int n, double * d, double ** v, int * nrot);


/* Gauss-Jordan elimination.
 *
 * To invert the matrix m[][], destroying m[][] and leaving the inverse
 * in its place, you can do
 *           gauss_jordan(m, n, (double **)0, 0);
 *
 * To see how the other arguments are used in solving linear systems,
 * see the Numerical Recipes book.
 */
extern void gauss_jordan(double ** m, int n, double ** b, int k);


/* LU Decomposition. Crout algorithm.
 *    m = L.U
 * The pointer arguments are assumed to be already allocated.
 * m is replaced by a matrix whose upper triangular part, including
 * the diagonal, is the matrix U, and whose off-diagonal lower
 * triangular part is the off-diagonal lower triangular part of
 * the matrix L. All the diagonal elements of L are equal to 1,
 * by construction, and do not appear in the result.
 */
extern void lu_decompose(double ** m, int n, int * indx, double * d);


/* LU Inversion.
 * Destroys the matrix m[][], placing the result in m_inv[][].
 */
extern void lu_invert(double ** m, double ** m_inverse, int n);


/* Invert using jacobi(). The input matrix is not altered.
 */
extern void jacobi_invert(const double ** m, double ** m_inverse, int n);


/*   Singular Value Decomposition routine. Based on the
 *   algorithm of J.C. Nash, Compact Numerical Methods for
 *   Computers (New York: Wiley and Sons, 1979), chapter 3.
 *
 *   Performs the singular value decomposition of real matrix A.
 *   thus finding orthogonal Q,R and diagonal S such that
 *             A.Q = R.S
 *
 *   ncol = number of columns of A
 *   nrow = number of rows of A
 *   S = diagonal ncol x ncol
 *   Q = orthogonal ncol x ncol
 *   R = orthogonal nrow x ncol
 *
 *   The algorithm destroys A, replacing it by R.
 *   Q is assumed to be previously allocated.
 */
extern void sv_decompose(double ** A, double ** Q, double * S, int nrow, int ncol);


/* Transpose a matrix in place. */
extern void transpose_d(double **, unsigned long size);
extern void transpose_i(int    **, unsigned long size);
extern void transpose_l(long   **, unsigned long size);
extern void transpose_f(float  **, unsigned long size);


/* Given eigenvalues d[] and associated eigenvectors in v[][],
 * sort them according to the values in d[].
 * The result is ordered lowest to highest, i.e. d[0] will
 * contain the lowest eigenvalue and v[i][0] will be the associated
 * eigenvector.
 *
 * The columns of v[][] are sorted, in accordance with
 * the convention of the result returned by jacobi().
 *
 * I have tested this on some 2x2's and 3x3's. (GJ)
 */
extern void eigsort(double * d, double ** v, int n);


/* Version of eigsort which sorts by absolute value.
 * The result is ordered lowest to highest, i.e. d[0] will
 * contain the eigenvalue of smallest absolute value,
 * and v[i][0] will be the associated eigenvector,
 *
 * I have tested this on some 2x2's and 3x3's. (GJ)
 */
extern void eigsort_abs(double * d, double ** v, int n);


/* Multiply matrices, putting result into a pre-allocated matrix. */
extern void matrix_multiply_i(int ** m_result,
			      int ** m1, int ** m2,
			      unsigned long size
			      );
extern void matrix_multiply_l(long ** m_result,
			      long ** m1, long ** m2,
			      unsigned long size
			      );
extern void matrix_multiply_f(float ** m_result,
			      float ** m1, float ** m2,
			      unsigned long size
			      );
extern void matrix_multiply_d(double ** m_result,
			      double ** m1, double ** m2,
			      unsigned long size
			      );


/* Find maximum and minimum elements of a matrix. */
extern void matrix_min_max(double ** m, unsigned long size, double * max, double * min);


/* Multiply a matrix by a scalar. */
extern void matrix_multiply_scalar(double ** m, unsigned long size, double s);


/* Normalize a matrix so that no element has absolute value greater than 1. */
extern void matrix_normalize(double ** m, unsigned long size);


/* Dump a matrix to a stream in raw form.
 */
extern void dump_matrix(FILE * fd, double ** m, int size, const char * format);


#endif /* MATRICES_H_ */
