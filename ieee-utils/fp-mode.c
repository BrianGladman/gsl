#include <config.h>

#if defined(HAVE_LINUX_IEEE_INTERFACE)
#include "fp-linux.c"
#elif defined(HAVE_HPUX_IEEE_INTERFACE)
#include "fp-hpux.c"
#else
#include "fp-unknown.c" 
#endif
