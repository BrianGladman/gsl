@BOTTOM@

/* Define if you have inline */
#undef HAVE_INLINE

/* Define if you have PACKAGE */
#undef PACKAGE

/* Define if you have VERSION */
#undef VERSION

/* Define if you have the ansi CLOCKS_PER_SEC clock rate */
#undef HAVE_CLOCKS_PER_SEC

/* Defined if configure has guessed a missing ansi CLOCKS_PER_SEC clock rate */
#undef HAVE_GUESSED_CLOCKS_PER_SEC

/* Use configure's best guess for CLOCKS_PER_SEC if it is unknown */
#ifndef HAVE_CLOCKS_PER_SEC
#define CLOCKS_PER_SEC HAVE_GUESSED_CLOCKS_PER_SEC
#endif

/* Defined if you have ansi EXIT_SUCCESS and EXIT_FAILURE in stdlib.h */
#undef HAVE_EXIT_SUCCESS_AND_FAILURE

/* Use 0 and 1 for EXIT_SUCCESS and EXIT_FAILURE if we don't have them */
#ifndef HAVE_EXIT_SUCCESS_AND_FAILURE
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif

/* Define one of these if you have a known IEEE arithmetic interface */
#undef HAVE_LINUX_IEEE_INTERFACE
#undef HAVE_SUNOS4_IEEE_INTERFACE
#undef HAVE_SOLARIS_IEEE_INTERFACE
#undef HAVE_HPUX_IEEE_INTERFACE
#undef HAVE_TRU64_IEEE_INTERFACE

/* Define this if printf can handle %Lf for long double */
#undef HAVE_PRINTF_LONGDOUBLE
