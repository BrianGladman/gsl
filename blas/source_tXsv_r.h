
  const int nonunit = (Diag == CblasNonUnit);
  size_t ix, jx;
  size_t i, j;

  if(N == 0) return;

  if(TransA == CblasNoTrans) {
    /* form  x := inv( A )*x */

    if(Uplo == CblasUpper) {
      /* backsubstitution */

      if(nonunit) {
        const size_t max_ix = incX * (N-1);
        X[max_ix] = X[max_ix]/ACCESS_UP(MATRIX_VAR_NAME,N,0,N-1,N-1);
      }

      ix = incX*(N-2);
      for(i=N-2; i>=0; i--) {
        BASE_TYPE tmp = X[ix];
	jx = ix + incX;
	for(j=i+1; j<GSL_MIN(N,i+KBAND+1); j++) {
	  tmp -= ACCESS_UP(MATRIX_VAR_NAME,N,0,i,j) * X[jx];
	  jx += incX;
	}
        X[ix] = tmp/ACCESS_UP(MATRIX_VAR_NAME,N,0,i,i);
	ix -= incX;
      }
    }
    else {
      /* forward substitution */

      if(nonunit) {
        X[0] = X[0]/ACCESS_LO(MATRIX_VAR_NAME,N,0,0,0);
      }

      ix = incX;
      for(i=1; i<N; i++) {
        BASE_TYPE tmp = X[ix];
	jx = 0;
	for(j=GSL_MAX(0,i-KBAND); j<i; j++) {
	  tmp -= ACCESS_LO(MATRIX_VAR_NAME,N,0,i,j) * X[jx];
	  jx += incX;
	}
        X[ix] = tmp/ACCESS_LO(MATRIX_VAR_NAME,N,0,i,i);
	ix += incX;
      }
    }
  }
  else {
    /* form  x := inv( A' )*x */

    if(Uplo == CblasUpper) {
      /* forward substitution */

      if(nonunit) {
        X[0] = X[0]/ACCESS_UP(MATRIX_VAR_NAME,N,0,0,0);
      }

      ix = incX;
      for(i=1; i<N; i++) {
        BASE_TYPE tmp = X[ix];
	jx = 0;
	for(j=GSL_MAX(0,i-KBAND); j<i; j++) {
	  tmp -= ACCESS_UP(MATRIX_VAR_NAME,N,0,j,i) * X[jx];
	  jx += incX;
	}
        X[ix] = tmp/ACCESS_UP(MATRIX_VAR_NAME,N,0,i,i);
	ix += incX;
      }
    }
    else {
      /* backsubstitution */

      if(nonunit) {
        const size_t max_ix = incX * (N-1);
        X[max_ix] = X[max_ix]/ACCESS_LO(MATRIX_VAR_NAME,N,0,N-1,N-1);
      }

      ix = incX*(N-2);
      for(i=N-2; i>=0; i--) {
        BASE_TYPE tmp = X[ix];
	jx = ix + incX;
	for(j=i+1; j<GSL_MIN(N,i+KBAND+1); j++) {
	  tmp -= ACCESS_LO(MATRIX_VAR_NAME,N,0,j,i) * X[jx];
	  jx += incX;
	}
        X[ix] = tmp/ACCESS_LO(MATRIX_VAR_NAME,N,0,i,i);
	ix -= incX;
      }
    }
  }
