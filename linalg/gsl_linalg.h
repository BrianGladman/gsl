/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_LINALG_H__
#define __GSL_LINALG_H__

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef enum {
  GSL_LA_MOD_NONE=0,
  GSL_LA_MOD_TRANSPOSE=1,
  GSL_LA_MOD_CONJUGATE=2
}
gsl_la_matrix_mod_t;


/* Simple implementation of matrix multiply.
 * Calculates C = A.B
 *
 * exceptions: GSL_EFAULT, GSL_EBADLEN
 */
int
gsl_la_matmult(const gsl_matrix * A, const gsl_matrix * B, gsl_matrix * C);


/* Simple implementation of matrix multiply.
 * Allows transposition of either matrix, so it
 * can compute A.B or Trans(A).B or A.Trans(B) or Trans(A).Trans(B)
 *
 * exceptions: GSL_EFAULT, GSL_EBADLEN
 */
int
gsl_la_matmult_mod(const gsl_matrix * A, gsl_la_matrix_mod_t modA,
                   const gsl_matrix * B, gsl_la_matrix_mod_t modB,
                   gsl_matrix * C);


/* Singular Value Decomposition
 *
 * exceptions: 
 */
int
gsl_la_decomp_SV_impl(gsl_matrix * A,
                      gsl_matrix * Q,
                      gsl_vector * S,
                      double tolerance);


/* LU Decomposition, Crout's method
 *
 * exceptions: 
 */
int
gsl_la_decomp_LU(gsl_matrix * matrix,
                 gsl_permutation * permutation,
                 int * signum);


/* Linear Solve Using LU Decomposition
 *
 * exceptions: 
 */
int
gsl_la_solve_LU_impl(const gsl_matrix     * lu_matrix,
                     const gsl_permutation * permutation,
                     const gsl_vector     * rhs,
                     gsl_vector           * solution);

int gsl_la_invert_LU_impl (const gsl_matrix     * lu_matrix,
                           const gsl_permutation * permutation,
                           gsl_matrix           * inverse);

double gsl_la_det_LU (gsl_matrix * lu_matrix, int signum);
double gsl_la_lndet_LU (gsl_matrix * lu_matrix);
int gsl_la_sgndet_LU (gsl_matrix * lu_matrix, int signum);

/* QR decomposition */

int
gsl_la_decomp_QR_impl(gsl_matrix * matrix, gsl_vector * rdiag);

int
gsl_la_solve_QR_impl(const gsl_matrix     * qr_matrix,
                     const gsl_vector     * rdiag,
                     const gsl_vector     * rhs,
                     gsl_vector           * solution);

int
gsl_la_qrsolve_QR_impl (const gsl_matrix * q, const gsl_matrix * r,
                        const gsl_vector * rhs,
                        gsl_vector * solution);

int
gsl_la_Rsolve_QR_impl (const gsl_matrix * qr_matrix,
                       const gsl_vector * rdiag,
                       gsl_vector * solution);

int 
gsl_la_update_QR_impl (gsl_matrix * q, gsl_matrix * r,
                       gsl_vector * w, const gsl_vector * v);

int
gsl_la_QTvec_QR_impl (const gsl_matrix * qr_matrix, 
                   const gsl_vector * v, 
                   gsl_vector * result);

int
gsl_la_unpack_QR_impl (const gsl_matrix * qr, const gsl_vector * rdiag,
                       gsl_matrix * q, gsl_matrix * r);


int
gsl_la_Rsolve_impl (const gsl_matrix * r, gsl_vector * solution);


/* Q R P^T decomposition */

int
gsl_la_decomp_QRPT_impl(gsl_matrix * matrix, 
                        gsl_vector * rdiag,
                        gsl_permutation * permutation, 
                        int * signum);

int
gsl_la_solve_QRPT_impl(const gsl_matrix     * qr_matrix,
                       const gsl_vector     * rdiag,
                       const gsl_permutation * permutation,
                       const gsl_vector     * rhs,
                       gsl_vector           * solution);

int
gsl_la_qrsolve_QRPT_impl (const gsl_matrix * q, const gsl_matrix * r,
                          const gsl_permutation * permutation,
                          const gsl_vector * rhs,
                          gsl_vector * solution);

int
gsl_la_Rsolve_QRPT_impl (const gsl_matrix * qr_matrix,
                       const gsl_vector * rdiag,
                       const gsl_permutation * permutation,
                       gsl_vector * solution);

int 
gsl_la_update_QRPT_impl (gsl_matrix * q, gsl_matrix * r,
                         const gsl_permutation * permutation,
                         gsl_vector * u, const gsl_vector * v);

/* Linear Solve Using Householder Transformations
 *
 * exceptions: 
 */
int
gsl_la_solve_HH_impl(gsl_matrix * matrix,
                     gsl_vector * vec);


/* Linear solve for a symmetric tridiagonal system.
 *
 * The input vectors represent the NxN matrix as follows:
 *
 *     diag[0]  offdiag[0]             0    ...
 *  offdiag[0]     diag[1]    offdiag[1]    ...
 *           0  offdiag[1]       diag[2]    ...
 *           0           0    offdiag[2]    ...
 *         ...         ...           ...    ...
 */
int
gsl_la_solve_symm_tridiag_impl(const gsl_vector * diag,
                               const gsl_vector * offdiag,
                               const gsl_vector * rhs,
                               gsl_vector * solution);


/* Linear solve for a symmetric cyclic tridiagonal system.
 *
 * The input vectors represent the NxN matrix as follows:
 *
 *      diag[0]  offdiag[0]             0   .....  offdiag[N-1]
 *   offdiag[0]     diag[1]    offdiag[1]   .....
 *            0  offdiag[1]       diag[2]   .....
 *            0           0    offdiag[2]   .....
 *          ...         ...
 * offdiag[N-1]         ...
 */
int
gsl_la_solve_symm_cyc_tridiag_impl(const gsl_vector * diag,
                                   const gsl_vector * offdiag,
                                   const gsl_vector * rhs,
                                   gsl_vector * solution);


#endif /* __GSL_LINALG_H__ */
