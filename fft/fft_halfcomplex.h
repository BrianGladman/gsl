int
  gsl_fft_halfcomplex_pass_2 (const double from[],
			      double to[],
			      const unsigned int product,
			      const unsigned int n,
			      const gsl_complex twiddle[]);

int
  gsl_fft_halfcomplex_pass_3 (const double from[], double to[],
			      const unsigned int product,
			      const unsigned int n,
			      const gsl_complex twiddle1[],
			      const gsl_complex twiddle2[]);

int
  gsl_fft_halfcomplex_pass_4 (const double from[],
			      double to[],
			      const unsigned int product,
			      const unsigned int n,
			      const gsl_complex twiddle1[],
			      const gsl_complex twiddle2[],
			      const gsl_complex twiddle3[]);

int
  gsl_fft_halfcomplex_pass_5 (const double from[],
			      double to[],
			      const unsigned int product,
			      const unsigned int n,
			      const gsl_complex twiddle1[],
			      const gsl_complex twiddle2[],
			      const gsl_complex twiddle3[],
			      const gsl_complex twiddle4[]);

int
  gsl_fft_halfcomplex_pass_6 (double *from, double *to,
			      unsigned int product, unsigned int n,
			      gsl_complex * twiddle1, gsl_complex * twiddle2,
			      gsl_complex * twiddle3, gsl_complex * twiddle4,
			      gsl_complex * twiddle5);

int
  gsl_fft_halfcomplex_pass_n (const double from[],
			      double to[],
			      const unsigned int factor,
			      const unsigned int product,
			      const unsigned int n,
			      const gsl_complex twiddle[]);


