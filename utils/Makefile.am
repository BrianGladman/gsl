noinst_LTLIBRARIES = libutils.la

# Don't need to list alloca.c, etc., Automake includes them.
libutils_la_SOURCES = system.h placeholder.c
libutils_la_LIBADD = @LIBOBJS@ 
libutils_la_DEPENDENCIES = $(libutils_la_LIBADD)

EXTRA_DIST = README
