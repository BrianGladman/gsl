#ifndef _GSL_COMPLEX_H
#define _GSL_COMPLEX_H

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


#define  GSL_COMPLEX_REAL(z)     (z.dat[0])
#define  GSL_COMPLEX_IMAG(z)     (z.dat[1])
#define  GSL_COMPLEX_P_REAL(zp)  (zp->dat[0])
#define  GSL_COMPLEX_P_IMAG(zp)  (zp->dat[1])


void   gsl_complex_set(gsl_complex * z, double re, double im);
double gsl_complex_real(const gsl_complex * z);
double gsl_complex_imag(const gsl_complex * z);

void   gsl_complex_float_set(gsl_complex_float * z, float re, float im);
double gsl_complex_float_real(const gsl_complex_float * z);
double gsl_complex_float_imag(const gsl_complex_float * z);


#ifdef HAVE_INLINE

inline
void
gsl_complex_set(gsl_complex * z, double re, double im)
{
  GSL_COMPLEX_P_REAL(z) = re;
  GSL_COMPLEX_P_IMAG(z) = im;
}

inline
double
gsl_complex_real(const gsl_complex * z)
{
  return GSL_COMPLEX_P_REAL(z);
}

inline
double
gsl_complex_imag(const gsl_complex * z)
{
  return GSL_COMPLEX_P_IMAG(z);
}

#endif /* HAVE_INLINE */


#endif /* _GSL_COMPLEX_H */
