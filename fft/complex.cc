#include <complex>
#include <vector>

extern "C" {
#include <gsl_fft_complex.h>
} ;

class FFTServer {

private:
  gsl_fft_complex_wavetable _wavetable;
  _construct(unsigned int n) ;
  _destruct() ;

public:
  FFTServer();                 // default constructor, no FFT length assigned
  FFTServer(unsigned int n);   // constructor for FFT of length n
  FFTServer(FFTServer &w);     // copy constructor, explicit mem management
  operator= (FFTServer &w);    // assignment, explicit mem management
  ~FFTServer() ;               // explicit memory management

  void set_order(unsigned int n) ;

  vector<complex<double> > forward_fft (vector<complex<double> > & v) ;
  vector<complex<double> > backward_fft (vector<complex<double> > & v) ;
  vector<complex<double> > inverse_fft (vector<complex<double> > & v) ;

  forward_fft_in_place (vector<complex<double> > & v) ;
  backward_fft_in_place (vector<complex<double> > & v) ;
  inverse_fft_in_place (vector<complex<double> > & v) ;

} ;

extern "C" {
#include <gsl_errno.h>
} ;

#include <stdexcept>

class gsl_exception: public exception { } ;
class gsl_logic_error: public gsl_exception { } ;

class gsl_bad_alloc: public bad_alloc, public gsl_exception { } ;

class gsl_out_of_range: public out_of_range, public gsl_logic_error {
public:
  gsl_out_of_range (const string& what_arg): out_of_range (what_arg) { }
} ;

class gsl_length_error: public length_error, public gsl_logic_error {
public:
  gsl_length_error (const string& what_arg): length_error (what_arg) { }
} ;

  
void gsl_check_exceptions(int status)
{
  switch(status) 
    {
    case 0:
      return ;
      break ;

    case GSL_ENOMEM:
      throw gsl_bad_alloc() ;
      break ;

    case GSL_ERANGE:
      throw gsl_out_of_range("") ;
      break ;
      
    default:
      throw gsl_exception() ;
    }
} ;

FFTServer::_construct(unsigned int n)

{  
  int status = gsl_fft_complex_wavetable_alloc (n, &_wavetable);

  gsl_check_exceptions(status) ;

  status = gsl_fft_complex_init (n, &_wavetable);

  gsl_check_exceptions(status) ;

  status = gsl_fft_complex_generate_wavetable (n, &_wavetable);

  gsl_check_exceptions(status) ;

} ;


FFTServer::_destruct()
{
  int status = gsl_fft_complex_wavetable_free (&_wavetable);

  gsl_check_exceptions(status) ;
} ;


FFTServer::FFTServer()
{
  _wavetable.n = 0 ; // do nothing, I guess
} ;


FFTServer::FFTServer(unsigned int n)
{
  _construct(n) ;
} ;

FFTServer::FFTServer(FFTServer &rhs)
{
  if (this == &rhs) 
    return ;

  int n = rhs._wavetable.n ;
  _wavetable = rhs._wavetable ;

  int status = gsl_fft_complex_wavetable_alloc (n, &_wavetable);

  gsl_check_exceptions(status) ;

  status = gsl_fft_complex_wavetable_cpy (&_wavetable, &rhs._wavetable) ;

  gsl_check_exceptions(status) ;
}

FFTServer::~FFTServer()
{
  _destruct() ;
}

void FFTServer::set_order(unsigned int n)
{
  if (_wavetable.n == n) 
    return ;  // size is ok

  _destruct() ;
  _construct(n) ;
}

vector<complex<double> > 
FFTServer::forward_fft(vector<complex<double> > & v)
{
  vector<complex<double> > w = v ;
    
  if (_wavetable.n != v.size())  // resize, or create, if necessary
    {
      set_order(v.size()) ; 
    } ;
  
  int status = gsl_fft_complex_forward ((gsl_complex *) w.begin(), 
					v.size(), 
					&_wavetable) ;

  gsl_check_exceptions(status) ;

  return w ;
}

vector<complex<double> > 
FFTServer::backward_fft(vector<complex<double> > & v)
{
  vector<complex<double> > w = v ;

  if (_wavetable.n != v.size())  // resize, or create, if necessary
    {
      set_order(v.size()) ; 
    } ;
  
  int status = gsl_fft_complex_backward((gsl_complex *) w.begin(), 
					v.size(), 
					&_wavetable);

  gsl_check_exceptions(status) ;

  return w ;
}

vector<complex<double> > 
FFTServer::inverse_fft(vector<complex<double> > & v)
{
  vector<complex<double> > w = v ;

  if (_wavetable.n != v.size())  // resize, or create, if necessary
    {
      set_order(v.size()) ; 
    } ;
  
  int status = gsl_fft_complex_inverse((gsl_complex *) w.begin(), 
				       v.size(), 
				       &_wavetable) ;
  gsl_check_exceptions(status) ;
  
  return w ;
}


FFTServer::forward_fft_in_place(vector<complex<double> > & v)
{
  if (_wavetable.n != v.size())  // resize, or create, if necessary
    {
      set_order(v.size()) ; 
    } ;
  
  int status = gsl_fft_complex_forward ((gsl_complex *) v.begin(), 
					v.size(), 
					&_wavetable) ;
  
  gsl_check_exceptions(status) ;
}

FFTServer::backward_fft_in_place(vector<complex<double> > & v)
{
  if (_wavetable.n != v.size())  // resize, or create, if necessary
    {
      set_order(v.size()) ; 
    } ;
  
  int status = gsl_fft_complex_backward ((gsl_complex *) v.begin(), 
					v.size(), 
					&_wavetable) ;
  
  gsl_check_exceptions(status) ;
}

FFTServer::inverse_fft_in_place(vector<complex<double> > & v)
{
  if (_wavetable.n != v.size())  // resize, or create, if necessary
    {
      set_order(v.size()) ; 
    } ;
  
  int status = gsl_fft_complex_inverse ((gsl_complex *) v.begin(), 
					v.size(), 
					&_wavetable) ;
  
  gsl_check_exceptions(status) ;
}
