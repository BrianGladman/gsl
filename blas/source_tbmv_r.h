
  int nounit = ( Diag == CblasNonUnit );
  size_t i;
  size_t j;

  if(TransA == CblasNoTrans) {
    /* form  x := A*x */

    if(Uplo == CblasUpper) {
      size_t ix = 0;
      for(i=0; i<N; i++) {
        BASE_TYPE temp = 0.0;
        size_t jx = (i+1)*incX;
        for(j=i+1; j<GSL_MIN(N,i+K+1); j++) {
	  temp += X[jx] * A[lda * i + j];
	  jx += incX;
	}
	if(nounit) {
	  X[ix] = temp + X[ix] * A[lda * i + i];
	}
	else {
	  X[ix] += temp;
	}
	ix += incX;
      }
    }
    else {
      size_t ix = 0;
      for(i=N-1; i+1>=1; i--) {
        BASE_TYPE temp = 0.0;
	const size_t j_min = ( K>i ? 0 : i-K );
	size_t jx = j_min * incX;
        for(j=j_min; j<i; j++) {
	  temp += X[jx] * A[lda * i + j];
	  jx += incX;
	}
	if(nounit) {
	  X[ix] = temp + X[ix] * A[lda * i + i];
	}
	else {
	  X[ix] += temp;
	}
	ix -= incX;
      }
    }
  }
  else {
    /* form  x := A'*x */
    /* FIXME: bad access */

    if(Uplo == CblasUpper) {
      size_t ix = 0;
      for(i=0; i<N; i++) {
        BASE_TYPE temp = 0.0;
        size_t jx = (i+1)*incX;
        for(j=i+1; j<GSL_MIN(N,i+K+1); j++) {
	  temp += X[jx] * A[lda * j + i];
	  jx += incX;
	}
	if(nounit) {
	  X[ix] = temp + X[ix] * A[lda * i + i];
	}
	else {
	  X[ix] += temp;
	}
	ix += incX;
      }
    }
    else {
      size_t ix = 0;
      for(i=N-1; i+1>=1; i--) {
        BASE_TYPE temp = 0.0;
	const size_t j_min = ( K>i ? 0 : i-K );
        for(j=j_min; j<i; j++) {
	  temp += X[j * incX] * A[lda * j + i];
	}
	if(nounit) {
	  X[ix] = temp + X[ix] * A[lda * i + i];
	}
	else {
	  X[ix] += temp;
	}
	ix += incX;
      }
    }
  }
