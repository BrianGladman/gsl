/* Some macros for accessing matrix elements
 * in raw array storage blocks. Oh the torture...
 */
#ifndef _MATRIX_ACCESS_H
#define _MATRIX_ACCESS_H

#include <complex_internal.h>


/* number of elements in rows up to and including
 * row i, in a packed triangular matrix
 */
#define PACKED_TR_COUNT(N,i)  ((((i)+1)*(2*(N)-(i)))/2)

/* access standard dense format */
#define M_STANDARD_INDEX(N, lda, i, j)       ((lda)*(i) + (j))
#define M_STANDARD_ACCESS(A, N, lda, i, j)   (A[M_STANDARD_INDEX(N, lda, i, j)])

/* access standard dense format, complex */
#define M_STANDARD_ACCESS_CR(A, N, lda, i, j)   REAL(A, 1, M_STANDARD_INDEX(N, lda, i, j))
#define M_STANDARD_ACCESS_CI(A, N, lda, i, j)   IMAG(A, 1, M_STANDARD_INDEX(N, lda, i, j))



/* access packed upper triangular format */
#define M_PACKEDTRUP_INDEX(N, lda, i, j)       (PACKED_TR_COUNT(N,(i)-1) + (j)-(i))
#define M_PACKEDTRUP_ACCESS(A, N, lda, i, j)   (A[M_PACKEDTRUP_INDEX(N, lda, i, j)])

/* access packed lower triangular format */
#define M_PACKEDTRLO_INDEX(N, lda, i, j)       (((i)*((i)+1))/2 + (j))
#define M_PACKEDTRLO_ACCESS(A, N, lda, i, j)   (A[M_PACKEDTRLO_INDEX(N, lda, i, j)])


/* access packed upper triangular format, complex */
#define M_PACKEDTRUP_ACCESS_CR(A, N, lda, i, j)   REAL(A, 1, M_PACKEDTRUP_INDEX(N, lda, i, j))
#define M_PACKEDTRUP_ACCESS_CI(A, N, lda, i, j)   IMAG(A, 1, M_PACKEDTRUP_INDEX(N, lda, i, j))


/* access packed lower triangular format, complex */
#define M_PACKEDTRLO_ACCESS_CR(A, N, lda, i, j)   REAL(A, 1, M_PACKEDTRLO_INDEX(N, lda, i, j))
#define M_PACKEDTRLO_ACCESS_CI(A, N, lda, i, j)   IMAG(A, 1, M_PACKEDTRLO_INDEX(N, lda, i, j))


#endif  /* !_MATRIX_ACCESS_H */
