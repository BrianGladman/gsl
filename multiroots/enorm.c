static double enorm (const gsl_vector * f);

static double enorm (const gsl_vector * f) {
  double e2 = 0 ;
  size_t i, n = f->size ;
  for (i = 0; i < n ; i++) {
    double fi= gsl_vector_get(f, i);
    e2 += fi * fi ;
  }
  return sqrt(e2);
}

static double scaled_enorm (const gsl_vector * d, const gsl_vector * f);

static double scaled_enorm (const gsl_vector * d, const gsl_vector * f) {
  double e2 = 0 ;
  size_t i, n = f->size ;
  for (i = 0; i < n ; i++) {
    double fi= gsl_vector_get(f, i);
    double di= gsl_vector_get(d, i);
    double u = di * fi;
    e2 += u * u ;
  }
  return sqrt(e2);
}


