pkglib_LIBRARIES = libgslhistogram.a 

pkginclude_HEADERS = gsl_histogram.h gsl_histogram2d.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../block -I$(srcdir)/../err

bin_PROGRAMS = gsl-histogram
check_PROGRAMS = test test2d test_resample test2d_resample test_trap test2d_trap

TESTS = test test2d test_resample test2d_resample test_trap test2d_trap test_gsl_histogram.sh

gsl_histogram_SOURCES = gsl-histogram.c
gsl_histogram_LDADD = libgslhistogram.a ../block/libgslblock.a  ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
noinst_HEADERS = urand.c

test_LDADD = libgslhistogram.a ../block/libgslblock.a  ../err/libgslerr.a ../err/libgsltest.a
test2d_LDADD = libgslhistogram.a ../block/libgslblock.a  ../err/libgslerr.a ../err/libgsltest.a
test_resample_LDADD = libgslhistogram.a ../block/libgslblock.a ../err/libgslerr.a ../err/libgsltest.a
test2d_resample_LDADD = libgslhistogram.a ../block/libgslblock.a  ../err/libgslerr.a ../err/libgsltest.a
test_trap_LDADD = libgslhistogram.a ../block/libgslblock.a ../err/libgslerr.a ../err/libgsltest.a
test2d_trap_LDADD = libgslhistogram.a ../block/libgslblock.a  ../err/libgslerr.a ../err/libgsltest.a

test_DEPENDENCIES = libgslhistogram.a
test2d_DEPENDENCIES = libgslhistogram.a
test_resample_DEPENDENCIES = libgslhistogram.a
test2d_resample_DEPENDENCIES = libgslhistogram.a
test_trap_DEPENDENCIES = libgslhistogram.a
test2d_trap_DEPENDENCIES = libgslhistogram.a

test_SOURCES = test.c
test2d_SOURCES = test2d.c
test_resample_SOURCES = test_resample.c
test2d_resample_SOURCES = test2d_resample.c
test_trap_SOURCES = test_trap.c
test2d_trap_SOURCES = test2d_trap.c

CLEANFILES = test.txt test.dat

libgslhistogram_a_SOURCES = \
	add.c find.c get.c init.c params.c reset.c file.c pdf.c \
	gsl_histogram.h \
	add2d.c find2d.c get2d.c init2d.c params2d.c reset2d.c \
	file2d.c pdf2d.c \
	gsl_histogram2d.h 




