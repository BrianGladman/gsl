#ifndef __GSL_FFT_COMPLEX_H__
#define __GSL_FFT_COMPLEX_H__

#include <stddef.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_fft.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/*  Power of 2 routines  */


int gsl_fft_complex_radix2_forward (gsl_complex_packed_array data,
                                    const size_t stride,
                                    const size_t n);

int gsl_fft_complex_radix2_backward (gsl_complex_packed_array data,
                                     const size_t stride,
                                     const size_t n);

int gsl_fft_complex_radix2_inverse (gsl_complex_packed_array data,
                                    const size_t stride,
                                    const size_t n);

int gsl_fft_complex_radix2_transform (gsl_complex_packed_array data,
                                      const size_t stride,
                                      const size_t n,
                                      const gsl_fft_direction sign);

int gsl_fft_complex_radix2_dif_forward (gsl_complex_packed_array data,
                                        const size_t stride,
                                        const size_t n);

int gsl_fft_complex_radix2_dif_backward (gsl_complex_packed_array data,
                                         const size_t stride,
                                         const size_t n);

int gsl_fft_complex_radix2_dif_inverse (gsl_complex_packed_array data,
                                        const size_t stride,
                                        const size_t n);

int gsl_fft_complex_radix2_dif_transform (gsl_complex_packed_array data,
                                          const size_t stride,
                                          const size_t n,
                                          const gsl_fft_direction sign);

int gsl_fft_complex_bitreverse_order (gsl_complex_packed_array data,
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
    double *scratch;
  }
gsl_fft_wavetable_complex;


gsl_fft_wavetable_complex *gsl_fft_complex_alloc (size_t n);

void gsl_fft_complex_free (gsl_fft_wavetable_complex * wavetable);

int gsl_fft_complex_memcpy (gsl_fft_wavetable_complex * dest,
                         gsl_fft_wavetable_complex * src);


int gsl_fft_complex_forward (gsl_complex_packed_array data,
                             const size_t stride,
                             const size_t n,
                             const gsl_fft_wavetable_complex * wavetable);

int gsl_fft_complex_backward (gsl_complex_packed_array data,
                              const size_t stride,
                              const size_t n,
                              const gsl_fft_wavetable_complex * wavetable);

int gsl_fft_complex_inverse (gsl_complex_packed_array data,
                             const size_t stride,
                             const size_t n,
                             const gsl_fft_wavetable_complex * wavetable);

int gsl_fft_complex_transform (gsl_complex_packed_array data,
                               const size_t stride, const size_t n,
                               const gsl_fft_wavetable_complex * wavetable,
                               const gsl_fft_direction sign);

__END_DECLS

#endif /* __GSL_FFT_COMPLEX_H__ */
