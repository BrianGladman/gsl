#ifndef _GSL_FFT_COMPLEX_H
#define _GSL_FFT_COMPLEX_H

#include <gsl_math.h>

#include <gsl_complex.h>
#include <gsl_fft.h>

typedef struct
  {
    unsigned int n;
    unsigned int nf;
    unsigned int factor[64];
    complex *twiddle[64];
    complex *trig;
    complex *scratch;
  }
gsl_fft_complex_wavetable;

  /*  Power of 2 routines  */

int
gsl_fft_complex_radix2_forward (complex data[],
				const unsigned int n);

int
gsl_fft_complex_radix2_backward (complex data[],
				 const unsigned int n);

int
gsl_fft_complex_radix2_inverse (complex data[],
				const unsigned int n);

int
gsl_fft_complex_radix2 (complex data[],
			const unsigned int n,
			const gsl_fft_direction sign);


int
gsl_fft_complex_radix2_dif_forward (complex data[],
				    const unsigned int n);

int
gsl_fft_complex_radix2_dif_backward (complex data[],
				     const unsigned int n);

int
gsl_fft_complex_radix2_dif_inverse (complex data[],
				    const unsigned int n);

int
gsl_fft_complex_radix2_dif (complex data[],
			    const unsigned int n,
			    const gsl_fft_direction sign);

int gsl_fft_binary_logn (const unsigned int n);

int gsl_fft_complex_bitreverse_order (complex data[], 
				      const unsigned int n,
				      const unsigned int n_bits);

  /*  Mixed Radix general-N routines  */

int gsl_fft_complex_forward (complex data[],
			     const unsigned int n,
			     const gsl_fft_complex_wavetable * wavetable);

int gsl_fft_complex_backward (complex data[],
			      const unsigned int n,
			      const gsl_fft_complex_wavetable * wavetable);

int gsl_fft_complex_inverse (complex data[],
			     const unsigned int n,
			     const gsl_fft_complex_wavetable * wavetable);

int
  gsl_fft_complex (complex data[],
		   const unsigned int n,
		   const gsl_fft_complex_wavetable * wavetable,
		   const gsl_fft_direction sign);

int
  gsl_fft_complex_init (const unsigned int n,
			gsl_fft_complex_wavetable * wavetable);

int
  gsl_fft_complex_generate_wavetable (unsigned int n,
				      gsl_fft_complex_wavetable * wavetable);

int
  gsl_fft_complex_wavetable_alloc (unsigned int n,
				   gsl_fft_complex_wavetable * wavetable);

int
  gsl_fft_complex_wavetable_free (gsl_fft_complex_wavetable * wavetable);

int
  gsl_fft_complex_pass_2 (const complex from[],
			  complex to[],
			  const gsl_fft_direction sign,
			  unsigned int product,
			  unsigned int n,
			  const complex twiddle[]);

int
  gsl_fft_complex_pass_3 (const complex from[],
			  complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int product,
			  const unsigned int n,
			  const complex twiddle1[],
			  const complex twiddle2[]);

int
  gsl_fft_complex_pass_4 (const complex from[],
			  complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int product,
			  const unsigned int n,
			  const complex twiddle1[],
			  const complex twiddle2[],
			  const complex twiddle3[]);

int
  gsl_fft_complex_pass_5 (const complex from[],
			  complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int product,
			  const unsigned int n,
			  const complex twiddle1[],
			  const complex twiddle2[],
			  const complex twiddle3[],
			  const complex twiddle4[]);

int
  gsl_fft_complex_pass_6 (const complex from[],
			  complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int product,
			  const unsigned int n,
			  const complex twiddle1[],
			  const complex twiddle2[],
			  const complex twiddle3[],
			  const complex twiddle4[],
			  const complex twiddle5[]);

int
  gsl_fft_complex_pass_n (complex from[],
			  complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int factor,
			  const unsigned int product,
			  const unsigned int n,
			  const complex twiddle[]);

#endif /* _GSL_FFT_COMPLEX_H */
