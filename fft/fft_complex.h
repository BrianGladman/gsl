int
  gsl_fft_complex_pass_2 (const gsl_complex from[],
			  gsl_complex to[],
			  const gsl_fft_direction sign,
			  unsigned int product,
			  unsigned int n,
			  const gsl_complex twiddle[]);

int
  gsl_fft_complex_pass_3 (const gsl_complex from[],
			  gsl_complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int product,
			  const unsigned int n,
			  const gsl_complex twiddle1[],
			  const gsl_complex twiddle2[]);

int
  gsl_fft_complex_pass_4 (const gsl_complex from[],
			  gsl_complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int product,
			  const unsigned int n,
			  const gsl_complex twiddle1[],
			  const gsl_complex twiddle2[],
			  const gsl_complex twiddle3[]);

int
  gsl_fft_complex_pass_5 (const gsl_complex from[],
			  gsl_complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int product,
			  const unsigned int n,
			  const gsl_complex twiddle1[],
			  const gsl_complex twiddle2[],
			  const gsl_complex twiddle3[],
			  const gsl_complex twiddle4[]);

int
  gsl_fft_complex_pass_6 (const gsl_complex from[],
			  gsl_complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int product,
			  const unsigned int n,
			  const gsl_complex twiddle1[],
			  const gsl_complex twiddle2[],
			  const gsl_complex twiddle3[],
			  const gsl_complex twiddle4[],
			  const gsl_complex twiddle5[]);

int
  gsl_fft_complex_pass_7 (const gsl_complex from[],
			  gsl_complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int product,
			  const unsigned int n,
			  const gsl_complex twiddle1[],
			  const gsl_complex twiddle2[],
			  const gsl_complex twiddle3[],
			  const gsl_complex twiddle4[],
			  const gsl_complex twiddle5[],
			  const gsl_complex twiddle6[]);


int
  gsl_fft_complex_pass_n (gsl_complex from[],
			  gsl_complex to[],
			  const gsl_fft_direction sign,
			  const unsigned int factor,
			  const unsigned int product,
			  const unsigned int n,
			  const gsl_complex twiddle[]);

