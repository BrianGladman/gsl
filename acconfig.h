@BOTTOM@

/* Define if you have the ansi CLOCKS_PER_SEC clock rate */
#undef HAVE_CLOCKS_PER_SEC

/* Defined if configure has guessed a missing ansi CLOCKS_PER_SEC clock rate */
#undef HAVE_GUESSED_CLOCKS_PER_SEC

/* Use configure's best guess for CLOCKS_PER_SEC if it is unknown */
#ifndef HAVE_CLOCKS_PER_SEC
#define CLOCKS_PER_SEC HAVE_GUESSED_CLOCKS_PER_SEC
#endif

/* Define if you have ansi rand() */
#undef HAVE_RAND

/* Define if you have ansi RAND_MAX */
#undef HAVE_RAND_MAX

/* Defined if configure has guessed a missing ansi RAND_MAX */
#undef HAVE_GUESSED_RAND_MAX

#ifndef HAVE_RAND_MAX
#define RAND_MAX HAVE_GUESSED_RAND_MAX
#endif
