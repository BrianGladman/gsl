#define OFFSET(N, incX) ((incX) > 0 ?  0 : ((N) - 1) * (-(incX)))
#define BLAS_ERROR(x)  abort(); 

#define CONJUGATE(x) ((x) == CblasConjTrans)
#define TRANSPOSE(x) ((x) == CblasTrans || (x) == CblasConjTrans)
#define UPPER(x) ((x) == CblasUpper)
#define LOWER(x) ((x) == CblasLower)

/* Handling of packed complex types... */

#define REAL(a,i) (((BASE *) a)[2*(i)])
#define IMAG(a,i) (((BASE *) a)[2*(i)+1])

#define REAL0(a) (((BASE *)a)[0])
#define IMAG0(a) (((BASE *)a)[1])

#define GB(KU,KL,lda,i,j) ((KU+1+(i-j))*lda + j)

#define TRCOUNT(N,i) ((((i)+1)*(2*(N)-(i)))/2)

#define TBUP(N,i,j) 
#define TBLO(N,i,j) 

#define TPUP(N,i,j) (TRCOUNT(N,(i)-1)+(j)-(i))
#define TPLO(N,i,j) (((i)*((i)+1))/2 + (j))
