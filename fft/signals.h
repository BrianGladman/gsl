int FUNCTION(fft_signal_complex,pulse) (size_t k, 
					size_t n,
					size_t stride,
					double z_real, double z_imag,
					BASE data[],
					BASE fft[]);

int FUNCTION(fft_signal_complex,constant) (size_t n,
					   size_t stride,
					   double z_real,
					   double z_imag,
					   BASE data[],
					   BASE fft[]);

int FUNCTION(fft_signal_complex,exp) (int k,
				      size_t n,
				      size_t stride,
				      double z_real,
				      double z_imag,
				      BASE data[],
				      BASE fft[]);


int FUNCTION(fft_signal_complex,exppair) (int k1,
					  int k2,
					  size_t n,
					  size_t stride,
					  double z1_real,
					  double z1_imag,
					  double z2_real,
					  double z2_imag,
					  BASE data[],
					  BASE fft[]);

int FUNCTION(fft_signal_complex,noise) (size_t n,
					size_t stride,
					BASE data[],
					BASE fft[]);

int FUNCTION(fft_signal_real,noise) (size_t n,
				     size_t stride,
				     BASE data[],
				     BASE fft[]);

