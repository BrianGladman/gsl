#ifndef _GSL_COMPLEX_H
#define _GSL_COMPLEX_H


/* two consecutive built-in types as a complex number */
typedef double *       gsl_complex_packed ;
typedef float *        gsl_complex_packed_float  ;
typedef long double *  gsl_complex_packed_long_double ;

/* 2N consecutive built-in types as N complex numbers */
typedef double *       gsl_complex_packed_array ;
typedef float *        gsl_complex_packed_array_float  ;
typedef long double *  gsl_complex_packed_array_long_double ;


typedef struct
  {
    long double dat[2];
  }
gsl_complex_long_double;

typedef struct
  {
    double dat[2];
  }
gsl_complex;

typedef struct
  {
    float dat[2];
  }
gsl_complex_float;

#define GSL_REAL(z)     ((z).dat[0])
#define GSL_IMAG(z)     ((z).dat[1])
#define GSL_COMPLEX_P_REAL(zp)  ((zp)->dat[0])
#define GSL_COMPLEX_P_IMAG(zp)  ((zp)->dat[1])
#define GSL_COMPLEX_EQ(z1,z2) ((z1).dat[0] == (z2).dat[0] && \
			       ((z1).dat[1] == (z2).dat[1]))
#define GSL_COMPLEX_SET(z,x,y) do {(z).dat[0]=(x); (z).dat[1]=(y);} while(0)

#endif /* _GSL_COMPLEX_H */
