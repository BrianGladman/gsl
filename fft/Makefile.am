lib_LIBRARIES = libgslfft.a 

noinst_PROGRAMS = test_trap test_radix2 test benchmark errs

include_HEADERS = gsl_fft.h gsl_fft_complex.h gsl_fft_halfcomplex.h gsl_fft_real.h gsl_fft_signals.h gsl_dft_complex.h

noinst_HEADERS = fft_complex.h fft_halfcomplex.h fft_real.h

INCLUDES = -I$(srcdir)/../err -I$(srcdir)/../utils

TESTS = test_trap test_radix2 test

test_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_trap_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
test_radix2_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
errs_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a
benchmark_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a ../utils/libutils.a

test_DEPENDENCIES = libgslfft.a
test_trap_DEPENDENCIES = libgslfft.a
test_radix2_DEPENDENCIES = libgslfft.a
benchmark_DEPENDENCIES = libgslfft.a


test_SOURCES = test.c
test_trap_SOURCES = test_trap.c
test_radix2_SOURCES = test_radix2.c
benchmark_SOURCES = benchmark.c 

libgslfft_a_SOURCES = \
compare.c compare.h \
dft_complex.c \
gsl_dft_complex.h \
 \
gsl_fft.h \
 \
complex.c \
gsl_fft_complex.h \
complex_init.c \
complex_pass_2.c \
complex_pass_3.c \
complex_pass_4.c \
complex_pass_5.c \
complex_pass_6.c \
complex_pass_7.c \
complex_pass_n.c \
complex_radix2.c \
 \
factorize.c \
factorize.h \
bitreverse.c \
bitreverse.h \
 \
halfcomplex.c \
gsl_fft_halfcomplex.h \
halfcomplex_init.c \
halfcomplex_pass_2.c \
halfcomplex_pass_3.c \
halfcomplex_pass_4.c \
halfcomplex_pass_5.c \
halfcomplex_pass_n.c \
halfcomplex_radix2.c \
halfcomplex_unpack.c \
 \
real.c \
gsl_fft_real.h \
real_init.c \
real_pass_2.c \
real_pass_3.c \
real_pass_4.c \
real_pass_5.c \
real_pass_n.c \
real_radix2.c \
real_unpack.c \
 \
signals.c

