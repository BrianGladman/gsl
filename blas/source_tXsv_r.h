/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */

  const int nonunit = (Diag == CblasNonUnit);
  size_t ix, jx;
  size_t i, j;
  size_t id;

  if(N == 0) return;

  if(TransA == CblasNoTrans) {
    /* form  x := inv( A )*x */

    if(Uplo == CblasUpper) {
      /* backsubstitution */

      if(nonunit) {
        const size_t max_ix = incX * (N-1);
        X[max_ix] = X[max_ix]/ACCESS_UP(MATRIX_VAR_NAME,N,LDA,N-1,N-1);
      }

      ix = incX*(N-2);
      for(id=0; id<N-1; id++) {
        BASE_TYPE tmp = X[ix];
        i = N-2-id;
	jx = ix + incX;
	for(j=i+1; j<GSL_MIN(N,i+KBAND+1); j++) {
	  const BASE_TYPE Aij = ACCESS_UP(MATRIX_VAR_NAME,N,LDA,i,j);
	  tmp -= Aij * X[jx];
	  jx += incX;
	}
	if(nonunit) {
          X[ix] = tmp/ACCESS_UP(MATRIX_VAR_NAME,N,LDA,i,i);
	}
	else {
	  X[ix] = tmp;
	}
	ix -= incX;
      }
    }
    else {
      /* forward substitution */

      if(nonunit) {
        X[0] = X[0]/ACCESS_LO(MATRIX_VAR_NAME,N,LDA,0,0);
      }

      ix = incX;
      for(i=1; i<N; i++) {
        BASE_TYPE tmp = X[ix];
	const size_t j0 = (i > KBAND ? i-KBAND : 0 );
	jx = 0;
	for(j=j0; j<i; j++) {
	  const BASE_TYPE Aij = ACCESS_LO(MATRIX_VAR_NAME,N,LDA,i,j);
	  tmp -= Aij * X[jx];
	  jx += incX;
	}
	if(nonunit) {
          X[ix] = tmp/ACCESS_LO(MATRIX_VAR_NAME,N,LDA,i,i);
	}
	else {
	  X[ix] = tmp;
	}
	ix += incX;
      }
    }
  }
  else {
    /* form  x := inv( A' )*x */

    if(Uplo == CblasUpper) {
      /* forward substitution */

      if(nonunit) {
        X[0] = X[0]/ACCESS_UP(MATRIX_VAR_NAME,N,LDA,0,0);
      }

      ix = incX;
      for(i=1; i<N; i++) {
        BASE_TYPE tmp = X[ix];
	const size_t j0 = ( i > KBAND ? i-KBAND : 0 );
	jx = 0;
	for(j=j0; j<i; j++) {
	  const BASE_TYPE Aji = ACCESS_UP(MATRIX_VAR_NAME,N,LDA,j,i);
	  tmp -= Aji * X[jx];
	  jx += incX;
	}
	if(nonunit) {
          X[ix] = tmp/ACCESS_UP(MATRIX_VAR_NAME,N,LDA,i,i);
	}
	else {
	  X[ix] = tmp;
	}
	ix += incX;
      }
    }
    else {
      /* backsubstitution */

      if(nonunit) {
        const size_t max_ix = incX * (N-1);
        X[max_ix] = X[max_ix]/ACCESS_LO(MATRIX_VAR_NAME,N,LDA,N-1,N-1);
      }

      ix = incX*(N-2);
      for(id=0; id<N-1; id++) {
        BASE_TYPE tmp = X[ix];
        i = N-2-id;
	jx = ix + incX;
	for(j=i+1; j<GSL_MIN(N,i+KBAND+1); j++) {
	  const BASE_TYPE Aji = ACCESS_LO(MATRIX_VAR_NAME,N,LDA,j,i);
	  tmp -= Aji * X[jx];
	  jx += incX;
	}
	if(nonunit) {
          X[ix] = tmp/ACCESS_LO(MATRIX_VAR_NAME,N,LDA,i,i);
	}
	else {
	  X[ix] = tmp;
	}
	ix -= incX;
      }
    }
  }
