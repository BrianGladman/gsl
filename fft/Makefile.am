noinst_LIBRARIES = libgslfft.a 

check_PROGRAMS = test_trap test_radix2 test benchmark errs

include_HEADERS = gsl_fft.h gsl_fft_complex.h gsl_fft_halfcomplex.h gsl_fft_real.h gsl_dft_complex.h

noinst_HEADERS = fft.h fft_complex.h fft_halfcomplex.h fft_real.h

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../vector -I$(srcdir)/../err

TESTS = test_trap test_radix2 test

test_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a
test_trap_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a
test_radix2_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a
errs_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a
benchmark_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a

test_DEPENDENCIES = libgslfft.a
test_trap_DEPENDENCIES = libgslfft.a
test_radix2_DEPENDENCIES = libgslfft.a
benchmark_DEPENDENCIES = libgslfft.a

test_SOURCES = test.c signals.c
test_trap_SOURCES = test_trap.c signals.c
test_radix2_SOURCES = test_radix2.c signals.c
benchmark_SOURCES = benchmark.c 

libgslfft_a_SOURCES = \
compare.c compare.h \
dft_complex.c \
gsl_dft_complex.h \
 \
gsl_fft.h \
 \
c.c \
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
real_unpack.c

