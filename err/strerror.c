#include <gsl_errno.h>

const char *
gsl_strerror (const int gsl_errno)
{
  switch(gsl_errno)
    {
    case GSL_EDOM:
      return "input domain error" ;
      break ;
    case GSL_ERANGE:
      return "output range error" ;
      break ;
    case GSL_EFAULT:
      return "invalid pointer" ;
      break ;
    case GSL_EINVAL:
      return "invalid argument supplied by user" ;
      break ;
    case GSL_EFAILED:
      return "generic failure" ;
      break ;
    case GSL_EFACTOR:
      return "factorization failed" ;
      break ;
    case GSL_ESANITY:
      return "sanity check failed - shouldn't happen" ;
      break ;
    case GSL_ENOMEM:
      return "malloc failed" ;
      break ;
    case GSL_EBADFUNC:
      return "problem with user-supplied function";
      break ;
    case GSL_ERUNAWAY:
      return "iterative process is out of control";
      break ;
    case GSL_ETIMEOUT:
      return "exceeded max number of iterations" ;
      break;
    case GSL_EZERODIV:
      return "tried to divide by zero" ;
      break;
    case GSL_ETOL:
      return "invalid tolerance or failed to reach tolerance" ;
      break;
    case GSL_EUNDRFLW:
      return "underflow" ;
      break ;
    case GSL_EOVRFLW:
      return "overflow" ;
      break;
    case GSL_ELOSS:
      return "loss of accuracy" ;
      break ;
    default:
      return "unknown error code" ;
    }
}
