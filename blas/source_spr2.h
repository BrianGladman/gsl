
  size_t i, j;

  if(Uplo == CblasUpper) {
    size_t k = 0;
    for(i=0; i<N; i++) {
      for(j=i; j<N; j++) {
        A[k] += alpha*(X[incX*i]*Y[incY*j] + X[incX*j]*Y[incY*i]);
	k++;
      }
    }
  }
  else {
    size_t k = 0;
    for(i=0; i<N; i++) {
      for(j=0; j<=i; j++) {
        A[k] += alpha*(X[incX*i]*Y[incY*j] + X[incX*j]*Y[incY*i]);
	k++;
      }
    }
  }
