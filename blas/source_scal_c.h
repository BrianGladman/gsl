
  size_t n;
  for(n=0; n<N; n++) {
    BASE_TYPE tmpR = REAL(X,incX,n);
    BASE_TYPE tmpI = IMAG(X,incX,n);
    REAL(X,incX,n) = tmpR*REAL0(alpha) - tmpI*IMAG0(alpha);
    IMAG(X,incX,n) = tmpR*IMAG0(alpha) + tmpI*REAL0(alpha);
  }
