#ifndef GSL_FFT_SIGNALS_H
#define GSL_FFT_SIGNALS_H

#include <gsl_math.h>

#include <gsl_complex.h>
#include <gsl_fft.h>

int gsl_fft_signal_complex_pulse (size_t k, size_t n,
				  double z_real, double z_imag,
				  double data[],
				  double fft[]);

int gsl_fft_signal_complex_constant (size_t n,
				     double z_real,
				     double z_imag,
				     double data[],
				     double fft[]);

int gsl_fft_signal_complex_exp (int k,
				size_t n,
				double z_real,
				double z_imag,
				double data[],
				double fft[]);


int gsl_fft_signal_complex_exppair (int k1,
				    int k2,
				    size_t n,
				    double z1_real,
				    double z1_imag,
				    double z2_real,
				    double z2_imag,
				    double data[],
				    double fft[]);

int gsl_fft_signal_complex_noise (size_t n,
				  double data[],
				  double fft[]);

int gsl_fft_signal_real_noise (size_t n,
			       double data[],
			       double fft[]);

#endif /* GSL_FFT_SIGNALS_H */
