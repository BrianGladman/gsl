noinst_LIBRARIES = libgslfft.a 

check_PROGRAMS = test #test_trap test_radix2 benchmark errs

include_HEADERS = gsl_fft.h gsl_fft_complex.h gsl_fft_halfcomplex.h gsl_fft_real.h gsl_dft_complex.h

noinst_HEADERS = c_pass.h hc_pass.h real_pass.h \
signals.h c_main.c c_init.c c_pass_2.c c_pass_3.c c_pass_4.c c_pass_5.c \
c_pass_6.c c_pass_7.c c_pass_n.c c_radix2.c bitreverse.c \
hc_init.c hc_pass_2.c hc_pass_3.c \
hc_pass_4.c hc_pass_5.c hc_pass_n.c hc_radix2.c hc_unpack.c \
real.c real_init.c real_pass_2.c real_pass_3.c \
real_pass_4.c real_pass_5.c real_pass_n.c real_radix2.c real_unpack.c \
compare_source.c dft_source.c hc_main.c real_main.c

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../vector -I$(srcdir)/../err

TESTS = test #test_radix2 test_trap

test_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a
#test_trap_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a
#test_radix2_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a
#errs_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a
#benchmark_LDADD = libgslfft.a ../err/libgslerr.a ../err/libgsltest.a

test_DEPENDENCIES = libgslfft.a
#test_trap_DEPENDENCIES = libgslfft.a
#test_radix2_DEPENDENCIES = libgslfft.a
#benchmark_DEPENDENCIES = libgslfft.a

test_SOURCES = test.c signals.c
#test_trap_SOURCES = test_trap.c signals.c
#test_radix2_SOURCES = test_radix2.c signals.c
#benchmark_SOURCES = benchmark.c 

libgslfft_a_SOURCES = compare.c compare.h dft.c \
gsl_dft_complex.h gsl_dft_complex_float.h \
gsl_fft.h \
gsl_fft_complex.h gsl_fft_complex_float.h \
gsl_fft_halfcomplex.h gsl_fft_halfcomplex_float.h \
gsl_fft_real.h gsl_fft_real_float.h \
c.c c_float.c  factorize.c factorize.h bitreverse.h \
hc.c real.c