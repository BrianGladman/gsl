#include <config.h>

#if defined(HAVE_SPARCLINUX_IEEE_INTERFACE)
#include "fp-sparclinux.c"
#elif defined(HAVE_M68KLINUX_IEEE_INTERFACE)
#include "fp-m68klinux.c"
#elif defined(HAVE_LINUX_IEEE_INTERFACE)
#include "fp-linux.c"
#elif defined(HAVE_HPUX_IEEE_INTERFACE)
#include "fp-hpux.c"
#elif defined(HAVE_SUNOS4_IEEE_INTERFACE)
#include "fp-sunos4.c"
#elif defined(HAVE_SOLARIS_IEEE_INTERFACE)
#include "fp-solaris.c"
#elif defined(HAVE_IRIX_IEEE_INTERFACE)
#include "fp-irix.c"
#elif defined(HAVE_AIX_IEEE_INTERFACE)
#include "fp-aix.c"
#elif defined(HAVE_TRU64_IEEE_INTERFACE)
#include "fp-tru64.c"
#else
#include "fp-unknown.c" 
#endif



