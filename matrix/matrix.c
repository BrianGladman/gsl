/* Created: [GJ] Sat Apr 27 02:04:31 EDT 1996
 */
#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include "memory.h"
#include "error.h"
#include "constants.h"
#include "sorting.h"
#include "matrices.h"


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void jacobi_diag(double ** a, int n, double * d, double ** v, int * nrot)
{
  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c;
  char buf[100];
  
  double * b = new_vector_d(n);
  double * z = new_vector_d(n);
  if( b == 0 || z == 0){
    char buff[100];
    sprintf(buff,"jacobi_diag: allocation failure");
    push_error(buff, Error_Alloc_);
    return;
  }

  for(ip=0; ip<n; ip++) {
    for(iq=0; iq<n; iq++) v[ip][iq]=0.0;
    v[ip][ip] = 1.0;
  }
  for(ip=0; ip<n; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++)
	sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      free_vector(z);
      free_vector(b);
      return;
    }
    if (i < 4)
      tresh = 0.2*sm/(n*n);
    else
      tresh = 0.0;
    for(ip=0; ip<n-1; ip++) {
      for(iq=ip+1; iq<n; iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
	    && fabs(d[iq])+g == fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if (fabs(h)+g == fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for(j=0; j<=ip-1; j++){
	    ROTATE(a,j,ip,j,iq)
	    }
	  for(j=ip+1; j<=iq-1; j++){
	    ROTATE(a,ip,j,j,iq)
	    }
	  for(j=iq+1; j<n; j++){
	    ROTATE(a,ip,j,iq,j)
	    }
	  for (j=0; j<n; j++){
	    ROTATE(v,j,ip,j,iq)
	    }
	  ++(*nrot);
	}
      }
    }
    for (ip=0; ip<n; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }

  sprintf(buf,"jacobi_diag: too many iterations");
  push_error(buf, Error_ConvFail_);
}
#undef ROTATE


#define TINY 1.0e-20;
void lu_decompose(double **a, int n, int *indx, double *d)
{
  int i, imax, j, k;
  double big, dum, sum, temp;
  double * vv = new_vector_d(n);
  
  *d = 1.0;
  for(i=0; i<n; i++) {
    big = 0.0;
    for(j=0; j<n; j++)
      if((temp=fabs(a[i][j])) > big) big = temp;
    if (big == 0.0) {
      char buff[100];
      sprintf(buff,"lu_decompose: singular matrix");
      push_error(buff, Error_ConvFail_);
      free_vector(vv);
      return;
    }
    vv[i] = 1.0/big;
  }
  for(j=0;j<n;j++) {
    for(i=0; i<j; i++) {
      sum = a[i][j];
      for(k=0; k<i; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for(i=j; i<n; i++){
      sum = a[i][j];
      for(k=0; k<j; k++)
	sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if((dum = vv[i]*fabs(sum)) >= big) {
	big = dum;
	imax = i;
      }
    }
    if(j != imax) {
      for(k=0; k<n; k++) {
	dum = a[imax][k];
	a[imax][k] = a[j][k];
	a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if(a[j][j] == 0.0) a[j][j] = TINY;
    if(j != n-1) {
      dum = 1.0/(a[j][j]);
      for(i=j+1; i<n; i++) a[i][j] *= dum;
    }
  }
  free_vector(vv);
}
#undef TINY


/* LU back substitution.
 * The arguments are assumed to be already allocated.
 */
static void lu_bksb(double **a, int n, int *indx, double b[])
{
  int i, ii=-1, ip, j;
  double sum;
  
  for(i=0; i<n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if(ii >= 0)
      for(j=ii; j<=i-1; j++) sum -= a[i][j]*b[j];
    else if(sum != 0.) ii=i;
    b[i] = sum;
  }
  for(i=n-1; i>=0; i--) {
    sum = b[i];
    for(j=i+1; j<n; j++) sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];
  }
}


void lu_invert(double **a, double **ainv, int n)
{
  int i, j;

  double d;
  double * col = new_vector_d(n);
  int * indx   = new_vector_i(n);

  lu_decompose(a,n,indx,&d);

  for(j=0; j<n; j++) {
    for(i=0; i<n; i++) { col[i] = 0.0; }
    col[j] = 1.0;
    
    lu_bksb(a,n,indx,col);

    for(i=0; i<n; i++) { ainv[i][j] = col[i]; }
  }

  free_vector(indx);
  free_vector(col);
}


#define SWAP(a,b) {temp=(a); (a)=(b); (b)=temp;}
void gauss_jordan(double **a, int n, double **b, int m)
{
  int i, icol, irow, j, k, l, ll;
  int * indxc = new_vector_i(n);
  int * indxr = new_vector_i(n);
  int * ipiv  = new_vector_i(n);
  double big, dum, pivinv, temp;

  for(j=0; j<n; j++) ipiv[j] = -1;

  for(i=0; i<n; i++){
    big = 0.;
    for(j=0; j<n; j++)
      if(ipiv[j] != 0)
	for(k=0; k<n; k++) {
	  if(ipiv[k] == -1) {
	    if(fabs(a[j][k]) >= big) {
	      big = fabs(a[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }
	  else if(ipiv[k] > 0) {
	    char buff[100];
	    sprintf(buff,"gaussjordan: singular matrix 0");
	    push_error(buff, Error_ConvFail_);
	    free_vector(indxc);
	    free_vector(indxr);
	    free_vector(ipiv);
	    return;
	  }
	}
    
    ++(ipiv[icol]);
    
    if(irow != icol){
      for(l=0; l<n; l++) SWAP(a[irow][l],a[icol][l])
      for(l=0; l<m; l++) SWAP(b[irow][l],b[icol][l])
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if(a[icol][icol] == 0.){
      char buff[100];
      sprintf(buff,"gaussjordan: singular matrix 1");
      push_error(buff, Error_ConvFail_);
      free_vector(indxc);
      free_vector(indxr);
      free_vector(ipiv);
      return;
    }
    pivinv = 1./a[icol][icol];
    a[icol][icol] = 1.;
    for(l=0; l<n; l++) a[icol][l] *= pivinv;
    for(l=0; l<m; l++) b[icol][l] *= pivinv;

    for(ll=0; ll<n; ll++)
      if(ll != icol){
	dum = a[ll][icol];
	a[ll][icol] = 0.;
	for(l=0; l<n; l++) a[ll][l] -= a[icol][l] * dum;
	for(l=0; l<m; l++) b[ll][l] -= b[icol][l] * dum;
      }
  }

  for(l=n-1; l>=0; l--) {
    if(indxr[l] != indxc[l])
      for(k=0; k<n; k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }

  free_vector(ipiv);
  free_vector(indxr);
  free_vector(indxc);
}
#undef SWAP


void jacobi_invert(const double ** a, double ** ainv, int n)
{
  int i,j,k,l;
  int nrot;

  double * eval = new_vector_d(n);
  double **evec = new_matrix_d(n);
  double **inv_diag = new_matrix_d(n);
  if(eval == (double *)0){
    char buff[100];
    sprintf(buff,"jacobi_invert: allocation failure");
    push_error(buff, Error_Alloc_);
    free_matrix_d(evec, n);
    free_matrix_d(inv_diag, n);
    return;
  }

  for(i=0; i<n; i++) {
    memcpy(ainv[i], a[i], n*sizeof(double));
  }

  jacobi_diag(ainv, n, eval, evec, &nrot);

  for(i=0; i<n; i++) {
    if(fabs(eval[i]) < 1.e-150){
      char buff[100];
      sprintf(buff,"jacobi_invert: apparent singularity");
      push_error(buff, Error_ConvFail_);
      free_vector(eval);
      free_matrix_d(evec, n);
      free_matrix_d(inv_diag, n);
      return;
    }
  }

  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      inv_diag[i][j] = ( i==j ? 1./eval[i] : 0.);
  
  for(i=0; i<n; i++)
    for(j=0; j<n; j++){
      ainv[i][j] = 0.;

      for(k=0; k<n; k++)
	for(l=0; l<n; l++)
	  ainv[i][j] += evec[i][l] * inv_diag[l][k] * evec[j][k];
    }

  free_matrix_d(inv_diag, n);
  free_matrix_d(evec, n);
  free_vector(eval);
}


#define GEN_TRANS(FUNC, TYPE)                     \
void FUNC(TYPE ** mat, unsigned long n)           \
{                                                 \
  int i,j;                                        \
  for(i=0; i<n; i++){                             \
    for(j=i+1; j<n; j++){                         \
      TYPE temp = mat[i][j];                      \
      mat[i][j] = mat[j][i];                      \
      mat[j][i] = temp;                           \
    }                                             \
  }                                               \
}                                                 \

GEN_TRANS(transpose_i, int)
GEN_TRANS(transpose_l, long)
GEN_TRANS(transpose_d, double)
GEN_TRANS(transpose_f, float)


void eigsort(double d[], double **v, int n)
{
  int k,j,i;
  double p;
  
  for(i=0; i<n-1; i++) {
    p=d[k=i];
    for(j=i+1; j<n; j++)
      if (d[j] <= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for(j=0; j<n; j++) {
        p=v[j][i];
        v[j][i]=v[j][k];
        v[j][k]=p;
      }
    }
  }
}


void eigsort_abs(double d[], double **v, int n)
{
  int k,j,i;
  double p;
  
  for(i=0;i<n-1;i++) {
    p=d[k=i];
    for(j=i+1; j<n; j++)
      if (fabs(d[j]) <= fabs(p)) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for(j=0; j<n; j++) {
        p=v[j][i];
        v[j][i]=v[j][k];
        v[j][k]=p;
      }
    }
  }
}

#define SVDTol_ 1.0e-8             /* tolerance value */
#define Macheps_ dbl_precision     /* machine precision */

void sv_decompose(double **A, double **Q, double *S, int nrow, int ncol)
{
  int i,j,k;

  /* Initialize the rotation counter and the sweep counter. */
  int count = 1;
  int sweep = 0;
  int sweepmax = ncol; 

  /* Always do at least 6 sweeps. */
  sweepmax = Max(sweepmax, 6);

  /* Set Q to the identity matrix. */
  for(i=0; i<ncol; i++){
    for(j=0; j<ncol; j++){
      Q[i][j] = 0.0;
    }
    Q[i][i] = 1.0;
  }

  /* Orthogonalize A by plane rotations. */
  while(count > 0 && sweep <= sweepmax){

    /* Initialize rotation counter. */
    count = ncol*(ncol-1)/2;        

    for(j=0; j<ncol-1; j++){
      for(k=j+1; k<ncol; k++){
        double p = 0.;
	double q = 0.;
	double r = 0.;
	double cosine, sine;
	double v;

        for(i=0;i<nrow;i++){
          p+=A[i][j]*A[i][k];    /* quantities in rotation angles */
          q+=A[i][j]*A[i][j];
          r+=A[i][k]*A[i][k];
        }
        if(q*r<Macheps_){count--; continue;}  /* column elements of A small */
        if(p*p/(q*r)<SVDTol_){count--; continue;} /* columns j,k orthogonal */ 
        if(q<r){cosine=0; sine=1;} else {
          q=q-r; 
          v=sqrt(4.0*p*p + q*q);
          cosine=sqrt((v+q)/(2.0*v));     /* rotation angles */
          sine=p/(v*cosine);
        }
        for(i=0;i<nrow;i++){     /* apply rotation to A */
          r=A[i][j];
          A[i][j]=r*cosine + A[i][k]*sine;
          A[i][k]= -r*sine + A[i][k]*cosine;
        }
        for(i=0;i<ncol;i++){     /* apply rotation to Q */
          r=Q[i][j];
          Q[i][j]=r*cosine + Q[i][k]*sine;
          Q[i][k]= -r*sine + Q[i][k]*cosine;
        }
      }
    }
    /* Sweep completed. */
    sweep++;
  }

  /* 
   * Orthogonalization complete. Compute singular values.
   */

  if(count > 0) {
    char buff[100];
    sprintf(buff,"sv_decompose: reached sweep limit"); 
    push_error(buff, Error_ConvFail_);
    return;
  }

  for(j=0;j<ncol;j++){
    double q = 0.0;

    /* Calculate singular values. */
    for(i=0; i<nrow; i++){
      q += A[i][j]*A[i][j];
    }
    S[j] = sqrt(q);

    /* Normalize vectors. */
    for(i=0; i<nrow; i++){
      A[i][j] /= S[j];
    }
  }
}
#undef SVDTol_
#undef Macheps_


/* There is only one thing we do to help the compiler here.
   The memory access to result[i][j] in the innermost step cannot
   generically be optimised away because of the possible pointer
   overlap problem. Therefore we move the access to a temp variable.
   This, together with compiler optimization to do some trivial
   code motion, provides near best possible performance.
 */
#define GEN_MATMULT(FUNC, TYPE)                                       \
void FUNC(TYPE ** result, TYPE ** m1, TYPE ** m2, unsigned long size) \
{                                                                     \
  if(m1 != 0 && m2 != 0) {                                            \
    unsigned long i, j, k;                                            \
    double temp;                                                      \
    for(i=0; i<size; i++){                                            \
      for(j=0; j<size; j++){                                          \
        temp = m1[i][0] * m2[0][j];				      \
	for(k=1; k<size; k++){                                        \
	  temp += m1[i][k] * m2[k][j];                                \
	}                                                             \
	result[i][j] = temp;                                          \
      }                                                               \
    }                                                                 \
  }								      \
}                                                                     \

GEN_MATMULT(matrix_multiply_i, int)
GEN_MATMULT(matrix_multiply_l, long)
GEN_MATMULT(matrix_multiply_f, float)
GEN_MATMULT(matrix_multiply_d, double)



void matrix_min_max(double ** m, unsigned long size, double * max, double * min)
{
  unsigned long i, j;
  double min_t = max_double;
  double max_t = -max_double;
  for(i=0; i<size; i++) {
    for(j=0; j<size; j++) {
      max_t = Max(max_t, m[i][j]);
      min_t = Min(min_t, m[i][j]);
    }
  }
  *min = min_t;
  *max = max_t;
}


void matrix_multiply_scalar(double ** m, unsigned long size, double s)
{
  unsigned long i, j;
  
  for(i=0; i<size; i++) {
    for(j=0; j<size; j++) {
      m[i][j] *= s;
    }
  }
}


void matrix_normalize(double ** m, unsigned long size)
{
  double max, min;
  double abs_max;
  matrix_min_max(m, size, &max, &min);
  abs_max = Max(fabs(max), fabs(min));
  matrix_multiply_scalar(m, size, 1./abs_max);
}


void dump_matrix(FILE * fd, double ** m, int size, const char * format)
{
  if(fd == 0) {
    push_error("dump_matrix: (nil) stream", Error_General_);
  }
  else {
    int i, j;
    for(i=0; i<size; i++) {
      for(j=0; j<size; j++) {
	fprintf(fd,format,m[i][j]);
	fprintf(fd," ");
      }
      fprintf(fd,"\n");
    }
  }
}
