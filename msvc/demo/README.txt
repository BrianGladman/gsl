Demonstration Workspace for GSL with Microsoft Visual C++
=========================================================

Workspace:  demo.dsw

Before compiling the project you need to check that the directory
where you installed GSL is present in the Microsoft Visual Studio
include and library paths.  The directory is added automatically to
these paths if you install the library with the automated setup
program.

The Visual C++ options that need to be present are,

    Tools --  Options -- Directories
          Include Files
             "C:\Program Files\GSL\include"
          Library Files: 
             "C:\Program Files\GSL\lib"

With these settings you should be able to compile and link the program
main.c in this directory.

You can choose from two project configurations, 'Release' and
'Debug'. The Release project configuration links with the libraries
libgsl.lib and libgslcblas.lib.  The Debug configuration links with
libgsld.lib and libgslcblasd.lib.

The output from the program should be the same as using GSL on Unix,

    bjg|debian> gcc main.c  -lgsl  -lgslcblas -lm
    bjg|debian> ./a.out 
    Here are ten random numbers in the range 0-99:
     66 36 72 68 57 81 27 83 13 95
    blas operation DAXPY
    x: 1 2 3 4 5  y: 5.5 4.4 3.3 2.2 1.1  a x + y: 6.5 6.4 6.3 6.2 6.1

If you see the same output, congratulations -- the library is working.

