lib_LTLIBRARIES = libgslfft.la 

check_PROGRAMS = test_trap test_radix2 test benchmark errs

include_HEADERS = gsl_fft.h gsl_fft_complex.h gsl_fft_halfcomplex.h gsl_fft_real.h gsl_fft_signals.h gsl_dft_complex.h

noinst_HEADERS = fft_complex.h fft_halfcomplex.h fft_real.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../err -I$(srcdir)/../utils

TESTS = test_trap test_radix2 test

test_LDADD = libgslfft.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.a
test_trap_LDADD = libgslfft.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.a
test_radix2_LDADD = libgslfft.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.a
errs_LDADD = libgslfft.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.a
benchmark_LDADD = libgslfft.la ../err/libgslerr.la ../err/libgsltest.la ../utils/libutils.a

test_DEPENDENCIES = libgslfft.la
test_trap_DEPENDENCIES = libgslfft.la
test_radix2_DEPENDENCIES = libgslfft.la
benchmark_DEPENDENCIES = libgslfft.la


test_SOURCES = test.c
test_trap_SOURCES = test_trap.c
test_radix2_SOURCES = test_radix2.c
benchmark_SOURCES = benchmark.c 

libgslfft_la_SOURCES = \
compare.c compare.h \
dft_complex.c \
gsl_dft_complex.h \
 \
gsl_fft.h \
 \
complex.c \
gsl_fft_complex.h \
c_init.c \
c_pass_2.c \
c_pass_3.c \
c_pass_4.c \
c_pass_5.c \
c_pass_6.c \
c_pass_7.c \
c_pass_n.c \
c_radix2.c \
 \
factorize.c \
factorize.h \
bitreverse.c \
bitreverse.h \
 \
hc.c \
gsl_fft_halfcomplex.h \
hc_init.c \
hc_pass_2.c \
hc_pass_3.c \
hc_pass_4.c \
hc_pass_5.c \
hc_pass_n.c \
hc_radix2.c \
hc_unpack.c \
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

