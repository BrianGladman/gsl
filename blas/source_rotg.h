
  BASE_TYPE roe   = ( fabs(*a) > fabs(*b) ? *a : *b );
  BASE_TYPE scale = fabs(*a) + fabs(*b);
  BASE_TYPE r, z;

  if( scale != 0.0 ) {
    BASE_TYPE aos = *a/scale;
    BASE_TYPE bos = *b/scale;
    r = scale * sqrt(aos*aos + bos*bos);
    r = sign(1.0,roe)*r;
    *c = *a/r;
    *s = *b/r;
    z = 1.0;
    if( fabs(*a) > fabs(*b) ) z = *s;
    if( fabs(*b) >= fabs(*a) && *c != 0.0 ) z = 1.0/(*c);
  }
  else {
   *c = 1.0;
   *s = 0.0;
   r = 0.0;
   z = 0.0;
  }
  
  *a = r;
  *b = z;
