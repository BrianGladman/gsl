#ifndef GSL_FFT_H
#define GSL_FFT_H

#include <gsl_complex.h>

#define REAL(a,stride,i) ((a)[2*(stride)*(i)])
#define IMAG(a,stride,i) ((a)[2*(stride)*(i)+1])

typedef enum
  {
    forward = -1, backward = +1   
  }
gsl_fft_direction;

/* this gives the sign in the formula

   h(f) = \sum x(t) exp(+/- 2 pi i f t) 
       
   where - is the forward transform direction and + the inverse direction */

#endif /* GSL_FFT_H */
