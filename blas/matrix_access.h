/* Some macros for accessing matrix elements
 * in raw array storage blocks. Oh the torture...
 */
#ifndef _MATRIX_ACCESS_H
#define _MATRIX_ACCESS_H


/* number of elements in rows up to and including
 * row i, in a packed triangular matrix
 */
#define PACKED_TR_COUNT(N,i)  ( (((i)+1)*(2*(N)-(i)))/2 )

/* access standard dense format */
#define M_STANDARD_ACCESS(A, N, lda, i, j)   (A[(lda)*(i) + (j)])

/* access packed triangular format */
#define M_PACKEDTRUP_ACCESS(A, N, lda, i, j)   (A[PACKED_TR_COUNT(N,(i)-1) + (j)-(i)])

/* access packed triangular format */
#define M_PACKEDTRLO_ACCESS(A, N, lda, i, j)   (A[((i)*((i)+1))/2 + (j)])


#endif  /* !_MATRIX_ACCESS_H */
