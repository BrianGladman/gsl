lib_LIBRARIES = libgslhistogram.a 

include_HEADERS = gsl_histogram.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../vector -I$(srcdir)/../err -I$(srcdir)/../utils

check_PROGRAMS = test test2d test_trap

TESTS = test test2d test_trap

test_LDADD = libgslhistogram.a ../vector/libgslvector.a  ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test2d_LDADD = libgslhistogram.a ../vector/libgslvector.a  ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_trap_LDADD = libgslhistogram.a ../vector/libgslvector.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a

test_DEPENDENCIES = libgslhistogram.a
test2d_DEPENDENCIES = libgslhistogram.a
test_trap_DEPENDENCIES = libgslhistogram.a

test_SOURCES = test.c
test2d_SOURCES = test2d.c
test_trap_SOURCES = test_trap.c

libgslhistogram_a_SOURCES = \
	add.c find.c get.c init.c params.c reset.c file.c pdf.c \
	gsl_histogram.h \
	init2d.c add2d.c file2d.c get2d.c params2d.c reset2d.c gsl_histgram2d.h 




