
  size_t i, j;
  const int nounit = ( Diag == CblasNonUnit );

  if(TransA == CblasNoTrans) {
    /* form  x:= A*x */

    if(Uplo == CblasUpper) {
      size_t k = 0;
      for(i=0; i<N; i++) {
        BASE_TYPE temp = ( nounit ? X[i*incX] * Ap[k] : X[i*incX] );
	k++;
        for(j=i+1; j<N; j++) {
          temp += Ap[k] * X[j * incX];
          k++;
        }
        X[i*incX] = temp;
      }
    }
    else {
      size_t k = (N*(N+1))/2 - 1;
      for(i=N-1; i>=0; i--) {
        BASE_TYPE temp = ( nounit ? X[i*incX] * Ap[k] : X[i*incX] );
	k--;
	for(j=GSL_MAX(0,i-1); j>=0; j--) {
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
      for(i=N-1; i>=0; i--) {
        size_t k = (i*(i+3))/2;
        BASE_TYPE temp = ( nounit ? X[i*incX] * Ap[k] : X[i*incX] );
	k -= N-i;
	for(j=GSL_MAX(0,i-1); j>=0; j--) {
          temp += Ap[k] * X[j * incX];
          k -= N-j;
        }
	X[i*incX] = temp;
      }
      
    }
    else {
      for(i=0; i<N; i++) {
        size_t k = (i*(i+3))/2;
        BASE_TYPE temp = ( nounit ? X[i*incX] * Ap[k] : X[i*incX] );
	k += i+1;
        for(j=i+1; j<N; j++) {
          temp += Ap[k] * X[j * incX];
          k += j+1;
        }
        X[i*incX] = temp;
      }
    }
  }
