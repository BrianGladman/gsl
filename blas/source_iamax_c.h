
  BASE_TYPE max = 0.0;
  CBLAS_INDEX n;
  CBLAS_INDEX i;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    const BASE_TYPE a = fabs(REAL(X, incX, i)) + fabs(IMAG(X, incX, i));
    if(a > max) {
      max = a;
      result = i;
    }
    i += incX;
  }
  return result;
