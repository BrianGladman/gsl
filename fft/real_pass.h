static void FUNCTION(fft_real,pass_2) (const BASE in[],
                                       const size_t istride,
                                       BASE out[],
                                       const size_t ostride,
                                       const size_t product,
                                       const size_t n,
                                       const gsl_complex twiddle[]);

static void FUNCTION(fft_real,pass_3) (const BASE in[], 
                                       const size_t istride,
                                       BASE out[],
                                       const size_t ostride,
                                       const size_t product,
                                       const size_t n,
                                       const gsl_complex twiddle1[],
                                       const gsl_complex twiddle2[]);

static void FUNCTION(fft_real,pass_4) (const BASE in[],
                                       const size_t istride,
                                       BASE out[],
                                       const size_t ostride,
                                       const size_t product,
                                       const size_t n,
                                       const gsl_complex twiddle1[],
                                       const gsl_complex twiddle2[],
                                       const gsl_complex twiddle3[]);

static void FUNCTION(fft_real,pass_5) (const BASE in[],
                                       const size_t istride,
                                       BASE out[],
                                       const size_t ostride,
                                       const size_t product,
                                       const size_t n,
                                       const gsl_complex twiddle1[],
                                       const gsl_complex twiddle2[],
                                       const gsl_complex twiddle3[],
                                       const gsl_complex twiddle4[]);

static void FUNCTION(fft_real,pass_n) (const BASE in[], 
                                       const size_t istride,
                                       BASE out[],
                                       const size_t ostride,
                                       const size_t factor,
                                       const size_t product,
                                       const size_t n,
                                       const gsl_complex twiddle[]);
