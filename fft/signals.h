int FUNCTION(fft_signal,complex_pulse) (size_t k, 
					size_t n,
					size_t stride,
					double z_real, double z_imag,
					BASE data[],
					BASE fft[]);

int FUNCTION(fft_signal,complex_constant) (size_t n,
					   size_t stride,
					   double z_real,
					   double z_imag,
					   BASE data[],
					   BASE fft[]);

int FUNCTION(fft_signal,complex_exp) (int k,
				      size_t n,
				      size_t stride,
				      double z_real,
				      double z_imag,
				      BASE data[],
				      BASE fft[]);


int FUNCTION(fft_signal,complex_exppair) (int k1,
					  int k2,
					  size_t n,
					  size_t stride,
					  double z1_real,
					  double z1_imag,
					  double z2_real,
					  double z2_imag,
					  BASE data[],
					  BASE fft[]);

int FUNCTION(fft_signal,complex_noise) (size_t n,
					size_t stride,
					BASE data[],
					BASE fft[]);

int FUNCTION(fft_signal,real_noise) (size_t n,
				     size_t stride,
				     BASE data[],
				     BASE fft[]);

