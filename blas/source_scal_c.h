
  size_t n;
  size_t ix = 0;
  for(n=0; n<N; n++) {
    BASE_TYPE tmpR = REAL(X,incX,ix);
    BASE_TYPE tmpI = IMAG(X,incX,ix);
    REAL(X,incX,ix) = tmpR*REAL0(alpha) - tmpI*IMAG0(alpha);
    IMAG(X,incX,ix) = tmpR*IMAG0(alpha) + tmpI*REAL0(alpha);
    ix += incX;
  }
