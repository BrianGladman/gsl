#include <config.h>

#if defined(HAVE_LINUX_IEEE_INTERFACE)
#include "fp-linux.c"
#elif defined(HAVE_HPUX_IEEE_INTERFACE)
#include "fp-hpux.c"
#elif defined(HAVE_SUNOS4_IEEE_INTERFACE)
#include "fp-sunos4.c"
#elif defined(HAVE_SOLARIS_IEEE_INTERFACE)
#include "fp-solaris.c"
#else
#include "fp-unknown.c" 
#endif
