
  size_t i, j;
  size_t ix, jx;
  const int nounit = ( Diag == CblasNonUnit );

  if(TransA == CblasNoTrans) {
    /* form  x:= A*x */

    if(Uplo == CblasUpper) {
      size_t k = 0;
      ix = 0;
      for(i=0; i<N; i++) {
        BASE_TYPE temp = ( nounit ? X[ix] * Ap[k] : X[ix] );
	k++;
	jx = (i+1)*incX;
        for(j=i+1; j<N; j++) {
          temp += Ap[k] * X[jx];
          k++;
	  jx += incX;
        }
        X[ix] = temp;
	ix += incX;
      }
    }
    else {
      size_t k = (N*(N+1))/2 - 1;
      for(i=N-1; i+1>=1; i--) {
        BASE_TYPE temp = ( nounit ? X[i*incX] * Ap[k] : X[i*incX] );
	const size_t j_min = ( i>1 ? i-1 : 0 );
	k--;
	for(j=j_min; j+1>=1; j--) {
          temp += Ap[k] * X[j * incX];
          k--;
        }
	X[i*incX] = temp;
      }
    }

  }
  else {
    /* form  x := A'*x */

    if(Uplo == CblasUpper) {
      for(i=N-1; i+1>=1; i--) {
        size_t k = (i*(i+3))/2;
        BASE_TYPE temp = ( nounit ? X[i*incX] * Ap[k] : X[i*incX] );
	const size_t j_min = ( i>1 ? i-1 : 0 );
	k -= N-i;
	for(j=j_min; j+1>=1; j--) {
          temp += Ap[k] * X[j * incX];
          k -= N-j;
        }
	X[i*incX] = temp;
      }
    }
    else {
      ix = 0;
      for(i=0; i<N; i++) {
        size_t k = (i*(i+3))/2;
        BASE_TYPE temp = ( nounit ? X[ix] * Ap[k] : X[ix] );
	k += i+1;
	jx = (i+1)*incX;
        for(j=i+1; j<N; j++) {
          temp += Ap[k] * X[jx];
          k += j+1;
	  jx += incX;
        }
        X[ix] = temp;
	ix += incX;
      }
    }
  }
