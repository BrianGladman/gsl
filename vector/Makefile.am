lib_LIBRARIES = libgslvector.a 

noinst_PROGRAMS = test

include_HEADERS = gsl_vector.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err -I$(srcdir)/../utils

TESTS = test

test_LDADD = libgslvector.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a

test_DEPENDENCIES = libgslvector.a

test_SOURCES = test.c

libgslvector_a_SOURCES = \
vector_init.c \
vector.c \
gsl_vector.h


