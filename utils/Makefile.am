pkglib_LIBRARIES = libutils.a

# Don't need to list alloca.c, etc., Automake includes them.
libutils_a_SOURCES = system.h placeholder.c
libutils_a_LIBADD = @LIBOBJS@ 
libutils_a_DEPENDENCIES = $(libutils_a_LIBADD)

EXTRA_DIST = README
