Compiling GSL with Microsoft Visual C++ 6.0
===========================================

1.  Build the workspace GSL.dsw.  This should produce the following
libraries,

        msvc/Debug/libgsl.lib
        msvc/Debug/libgslcblas.lib
or
        msvc/Release/libgsl.lib
        msvc/Release/libgslcblas.lib

depending on whether the Debug or Release build option is chosen.

[The header file msvc/config.h contains the appropriate configuration
for Microsoft Visual C++ and overrides the normal config.h.]

2.  Build the test workspace GSLTESTS.dsw

This should create the test programs in msvc/bin/*.exe.  To run the
tests use the batch file,

        MAKE_CHECK.bat

This produces an output log file "results.dat".  Any lines which don't
begin with PASS: indicate a problem.

3.  If the tests are successful install the files libgsl.lib,
libgslcblas.lib and gsl/gsl*.h.  You will need to add these to the
project settings for programs that you compile with GSL.

The library is built with the /LD or /LDd option.  This is compatible
with the default link option for single-threaded applications, /ML.
To use GSL in a multi-threaded application you may need to recompile
the library with another option, such as /MD.  See the Microsoft
Visual C++ Manual for details on link options.

The initial scripts were provided by José Miguel Buenaposada
(jmbuena@dia.fi.upm.es) and subsequently modified by Brian Gough
