#ifndef GSL_FFT_SIGNALS_H
#define GSL_FFT_SIGNALS_H

#include <gsl_math.h>

#include <gsl_complex.h>
#include <gsl_fft.h>

int
  gsl_fft_signal_complex_pulse (const size_t k,
				     const size_t n,
				     const double z_real,
				     const double z_imag,
				     gsl_complex data[],
				     gsl_complex fft[]);

int
  gsl_fft_signal_complex_constant (const size_t n,
					const double z_real,
					const double z_imag,
					gsl_complex data[],
					gsl_complex fft[]);

int
  gsl_fft_signal_complex_exp (const int k,
				   const size_t n,
				   const double z_real,
				   const double z_imag,
				   gsl_complex data[],
				   gsl_complex fft[]);


int
  gsl_fft_signal_complex_exppair (const int k1,
				       const int k2,
				       const size_t n,
				       const double z1_real,
				       const double z1_imag,
				       const double z2_real,
				       const double z2_imag,
				       gsl_complex data[],
				       gsl_complex fft[]);

int
  gsl_fft_signal_complex_noise (const size_t n,
				     gsl_complex data[],
				     gsl_complex fft[]);

int
  gsl_fft_signal_real_noise (const size_t n,
				  gsl_complex data[],
				  gsl_complex fft[]);

#endif /* GSL_FFT_SIGNALS_H */
