
  BASE_TYPE max = 0.0;
  CBLAS_INDEX n;
  CBLAS_INDEX i;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    if(fabs(X[i]) > max) {
      max = fabs(X[i]);
      result = i;
    }
    i += incX;
  }
  return result;
