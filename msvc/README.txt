This is a first attempt to make GSL compile on Visual C++ 6.0.

Initial scripts by José Miguel Buenaposada (jmbuena@dia.fi.upm.es)
Modified by Brian Gough

1. Start in the directory gsl/msvc/ and execute the batch file
./copy_gsl_headers.bat which copies the header files into the ../gsl/
subdirectory.

2.  Build the workspace GSL.dsw

The header file msvc/config.h contains the appropriate configuration
for Microsoft Visual C++ and overrides the normal config.h.  The same
is true for gsl_version.h which is not generated and the VC++ version
is in msvc/gsl_version.h.

This should produce the following libraries,

        msvc/Debug/libgsl.lib
        msvc/Debug/libgslcblas.lib
or
        msvc/Release/libgsl.lib
        msvc/Release/libgslcblas.lib

depending on whether the Debug or Release build option is chosen.

3.  Build the test workspace GSLTESTS.dsw

This should create the test programs in msvc/bin/*.exe.  To run the
tests use the batch file,

        MAKE_CHECK.bat

This produces an output log file "results.dat".  Any lines which don't
begin with PASS: indicate a problem.

4.  If the tests are successful install the files libgsl.lib,
libgslcblas.lib and ../gsl/gsl*.h
