#include <complex>
#include <iostream>
#include <vector>

extern "C" {
#include <gsl_fft_complex.h>
#include "bitreverse.h"
} ;


main() {
  double a[] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 } ;
  double b[] = { 0,0,0,0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0 } ;
  vector<complex<double> > x(16) ;

  for (int i = 0 ; i != x.size() ; i++) {
    *(x.begin() + i) = complex<double>(a[i],b[i]) ; 
    cout << *(x.begin() + i) << endl ;
  }

  gsl_fft_complex_bitreverse_order((gsl_complex *)x.begin(), 16, 4) ;

  for (int i = 0 ; i != x.size() ; i++) {
    cout << *(x.begin() + i) << endl ;
  }
}

class FFT {
private:
  gsl_fft_complex_wavetable _w;
public:
  FFT()                  // default constructor
  FFT(unsigned int n);   // constructor for length n
  FFT(FFT &w);           // copy constructor, explicit mem management
  operator= (FFT &w);    // assignment, explicit mem management
  ~FFT() ;               // explicit memory management
  forward() ;
  inverse() ;
  
}
  
FFT::FFT()
{
  // do nothing, I guess
}

FFT::FFT(unsigned int n)
{
  status = gsl_fft_complex_wavetable_alloc (n, &_w);
  status = gsl_fft_complex_init (n, &_w);
  status = gsl_fft_complex_generate_wavetable (n, &_w);
}

FFT::FFT(FFT &rhs)
{
  if (this == rhs) 
    return *this ;

  int n = rhs._w.n ;

  // no easy optimisation here, since the twiddle array is a set of
  // pointers into trig. Doing a memcpy on the members won't work --
  // it would require some manipulation to fix it up.

  status = gsl_fft_complex_wavetable_alloc (n, &_w); 
  status = gsl_fft_complex_init (n, &_w);
  status = gsl_fft_complex_generate_wavetable (n, &_w);  
}

FFT::~FFT()
{
  status = gsl_fft_complex_wavetable_free (n, &_w);
}

FFTComplex(vector<complex> v, FFT w, const gsl_fft_direction sign)
{
  int status = gsl_fft_complex(v.begin(), v.size(), w, sign) ;
  
  if (status) {
    cout << "wahhh" << status << endl ;
  }
}

