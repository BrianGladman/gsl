
  BASE_TYPE scale = 0.0;
  BASE_TYPE ssq   = 1.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    BASE_TYPE axi = fabs(REAL(X,incX,i));
    BASE_TYPE ayi = fabs(IMAG(X,incX,i));
    if(scale < axi) {
      ssq   = 1.0 + ssq*(scale/axi)*(scale/axi);
      scale = axi;
    }
    else {
      ssq += (axi/scale)*(axi/scale);
    }
    if(scale < ayi) {
      ssq   = 1.0 + ssq*(scale/ayi)*(scale/ayi);
      scale = axi;
    }
    else {
      ssq += (ayi/scale)*(ayi/scale);
    }
    i += incX;
  }
  return scale * sqrt(ssq);
