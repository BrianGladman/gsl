#ifndef GSL_FFT_COMPLEX_FLOAT_H
#define GSL_FFT_COMPLEX_FLOAT_H

#include <stddef.h>

#include <gsl_math.h>
#include <gsl_complex.h>
#include <gsl_fft.h>

/*  Power of 2 routines  */


int gsl_fft_complex_float_radix2_forward (gsl_complex_packed_array_float data,
                                          size_t stride,
                                          size_t n);

int gsl_fft_complex_float_radix2_backward (gsl_complex_packed_array_float data,
                                           size_t stride,
                                           size_t n);

int gsl_fft_complex_float_radix2_inverse (gsl_complex_packed_array_float data,
                                          size_t stride,
                                          size_t n);

int gsl_fft_complex_float_radix2_transform (gsl_complex_packed_array_float data,
                                            size_t stride,
                                            size_t n,
                                            gsl_fft_direction sign);

int gsl_fft_complex_float_radix2_dif_forward (gsl_complex_packed_array_float data,
                                              size_t stride,
                                              size_t n);

int gsl_fft_complex_float_radix2_dif_backward (gsl_complex_packed_array_float data,
                                               size_t stride,
                                               size_t n);

int gsl_fft_complex_float_radix2_dif_inverse (gsl_complex_packed_array_float data,
                                              size_t stride,
                                              size_t n);

int gsl_fft_complex_float_radix2_dif_transform (gsl_complex_packed_array_float data,
                                                size_t stride,
                                                size_t n,
                                                gsl_fft_direction sign);

int gsl_fft_complex_float_bitreverse_order (gsl_complex_packed_array_float data,
                                            size_t stride,
                                            size_t n,
                                            size_t n_bits);

/*  Mixed Radix general-N routines  */

typedef struct
  {
    size_t n;
    size_t nf;
    size_t factor[64];
    gsl_complex *twiddle[64];
    gsl_complex *trig;
    float *scratch;
  }
gsl_fft_wavetable_complex_float;


gsl_fft_wavetable_complex_float *gsl_fft_complex_float_alloc (size_t n);

void gsl_fft_complex_float_free (gsl_fft_wavetable_complex_float * wavetable);

int gsl_fft_complex_float_cpy (gsl_fft_wavetable_complex_float * dest,
                               gsl_fft_wavetable_complex_float * src);


int gsl_fft_complex_float_forward (gsl_complex_packed_array_float data,
                                   size_t stride,
                                   size_t n,
                                   const gsl_fft_wavetable_complex_float * wavetable);

int gsl_fft_complex_float_backward (gsl_complex_packed_array_float data,
                                    size_t stride,
                                    size_t n,
                                    const gsl_fft_wavetable_complex_float * wavetable);

int gsl_fft_complex_float_inverse (gsl_complex_packed_array_float data,
                                   size_t stride,
                                   size_t n,
                                   const gsl_fft_wavetable_complex_float * wavetable);

int gsl_fft_complex_float_transform (gsl_complex_packed_array_float data,
                                     size_t stride, size_t n,
                                     const gsl_fft_wavetable_complex_float * wavetable,
                                     gsl_fft_direction sign);

#endif /* GSL_FFT_COMPLEX_FLOAT_H */





