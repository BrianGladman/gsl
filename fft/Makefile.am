noinst_LIBRARIES = libfft.a 

bin_PROGRAMS = test_trap test_radix2 test benchmark errs

include_HEADERS = gsl_fft.h gsl_fft_complex.h gsl_fft_halfcomplex.h gsl_fft_real.h gsl_fft_signals.h gsl_dft_complex.h gsl_complex.h gsl_math.h

INCLUDES = -I$(srcdir)/../err

TESTS = test_trap test_radix2 test

test_LDADD = libfft.a ../err/liberr.a
test_trap_LDADD = libfft.a ../err/liberr.a
test_radix2_LDADD = libfft.a ../err/liberr.a
errs_LDADD = libfft.a ../err/liberr.a
benchmark_LDADD = libfft.a ../err/liberr.a

test_DEPENDENCIES = libfft.a
test_trap_DEPENDENCIES = libfft.a
test_radix2_DEPENDENCIES = libfft.a
benchmark_DEPENDENCIES = libfft.a

test_SOURCES = test.c
test_trap_SOURCES = test_trap.c
test_radix2_SOURCES = test_radix2.c
benchmark_SOURCES = benchmark.c 


libfft_a_SOURCES = \
autotest.c autotest.h \
compare.c compare.h \
gsl_complex.h \
dft_complex.c \
gsl_dft_complex.h \
 \
gsl_fft.h \
benchmark.c \
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
bitreverse.c \
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
signals.c \
 \
gsl_math.h
