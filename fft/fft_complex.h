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
  gsl_fft_complex_pass_7 (const complex from[],
			  complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int product,
			  const unsigned int n,
			  const complex twiddle1[],
			  const complex twiddle2[],
			  const complex twiddle3[],
			  const complex twiddle4[],
			  const complex twiddle5[],
			  const complex twiddle6[]);


int
  gsl_fft_complex_pass_n (complex from[],
			  complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int factor,
			  const unsigned int product,
			  const unsigned int n,
			  const complex twiddle[]);

