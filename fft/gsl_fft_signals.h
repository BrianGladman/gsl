#ifndef _GSL_FFT_SIGNALS_H
#define _GSL_FFT_SIGNALS_H

#include <gsl_math.h>

#include <gsl_complex.h>
#include <gsl_fft.h>

int
  gsl_fft_signal_complex_pulse (const unsigned int k,
				     const unsigned int n,
				     const double z_real,
				     const double z_imag,
				     gsl_complex data[],
				     gsl_complex fft[]);

int
  gsl_fft_signal_complex_constant (const unsigned int n,
					const double z_real,
					const double z_imag,
					gsl_complex data[],
					gsl_complex fft[]);

int
  gsl_fft_signal_complex_exp (const int k,
				   const unsigned int n,
				   const double z_real,
				   const double z_imag,
				   gsl_complex data[],
				   gsl_complex fft[]);


int
  gsl_fft_signal_complex_exppair (const int k1,
				       const int k2,
				       const unsigned int n,
				       const double z1_real,
				       const double z1_imag,
				       const double z2_real,
				       const double z2_imag,
				       gsl_complex data[],
				       gsl_complex fft[]);

int
  gsl_fft_signal_complex_noise (const unsigned int n,
				     gsl_complex data[],
				     gsl_complex fft[]);

int
  gsl_fft_signal_real_noise (const unsigned int n,
				  gsl_complex data[],
				  gsl_complex fft[]);

#endif /* _GSL_FFT_SIGNALS_H */
