Demonstration Workspace for GSL with Microsoft Visual C++
=========================================================

Workspace:  demo.dsw

The workspace is set up to look for include files in ..\include
and library files in ..\lib.

The 'Release' project configuration links with the libraries
libgsl.lib and libgslcblas.lib.  The 'Debug' configuration links with
libgsld.lib and libgslcblasd.lib.

If you can compile and link the program main.c in this directory then
your library is correctly installed.

The output from the program should be the same as using GSL on Unix,

    bjg|debian> gcc main.c  -lgsl  -lgslcblas -lm
    bjg|debian> ./a.out 
    Here are ten random numbers in the range 0-99:
     66 36 72 68 57 81 27 83 13 95
    blas operation DAXPY
    x: 1 2 3 4 5  y: 5.5 4.4 3.3 2.2 1.1  a x + y: 6.5 6.4 6.3 6.2 6.1
