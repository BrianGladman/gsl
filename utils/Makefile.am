# we don't install these, so we don't use libtool, at least that's 
# how mjb understands it.
#
noinst_LIBRARIES = libutils.a

# Don't need to list alloca.c, etc., Automake includes them.
libutils_a_SOURCES = system.h xmalloc.c xstrdup.c
libutils_a_LIBADD = @LIBOBJS@ @ALLOCA@
libutils_a_DEPENDENCIES = $(libutils_a_LIBADD)

EXTRA_DIST = README
