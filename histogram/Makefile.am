lib_LIBRARIES = libgslhistogram.a 

include_HEADERS = gsl_histogram.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../vector -I$(srcdir)/../err -I$(srcdir)/../utils

check_PROGRAMS = test #test_float test_int test_trap

TESTS = test #test_float test_int test_trap

test_LDADD = libgslhistogram.a ../vector/libgslvector.a  ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
#test_trap_LDADD = libgslvector.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a

test_DEPENDENCIES = libgslhistogram.a
#test_trap_DEPENDENCIES = libgslvector.a

test_SOURCES = test.c
#test_trap_SOURCES = test_trap.c

libgslhistogram_a_SOURCES = add.c find.c get.c init.c params.c reset.c file.c \
	 gsl_histogram.h




