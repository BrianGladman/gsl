
int
  gsl_fft_real_pass_2 (const double from[],
		       double to[],
		       const unsigned int product,
		       const unsigned int n,
		       const complex twiddle[]);

int
  gsl_fft_real_pass_3 (const double from[], double to[],
		       const unsigned int product,
		       const unsigned int n,
		       const complex twiddle1[],
		       const complex twiddle2[]);

int
  gsl_fft_real_pass_4 (const double from[],
		       double to[],
		       const unsigned int product,
		       const unsigned int n,
		       const complex twiddle1[],
		       const complex twiddle2[],
		       const complex twiddle3[]);

int
  gsl_fft_real_pass_5 (const double from[],
		       double to[],
		       const unsigned int product,
		       const unsigned int n,
		       const complex twiddle1[],
		       const complex twiddle2[],
		       const complex twiddle3[],
		       const complex twiddle4[]);

int
  gsl_fft_real_pass_6 (double *from, double *to,
		       unsigned int product, unsigned int n,
		       complex * twiddle1, complex * twiddle2,
		       complex * twiddle3, complex * twiddle4,
		       complex * twiddle5);

int
  gsl_fft_real_pass_n (const double from[], double to[],
		       const unsigned int factor,
		       const unsigned int product,
		       const unsigned int n,
		       const complex twiddle[]);
