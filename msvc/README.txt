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

For example, if the files are installed in the following locations,

       C:\gsl\include\gsl  - header files
       C:\gsl\lib          - lib files

the corresponding project settings should be,

       Project Settings
         C/C++
           Category: Preprocessor  
             Additional Include Directories: c:\gsl\include

       Project Settings
         Link
           Category: Input
             Object/Library Modules: libgsl.lib libgslcblas.lib ...
            and
             Additional Library Path: C:\gsl\lib

Make sure that the settings are applied to the appropriate
configuration (either 'Release' or 'Debug', depending on which version
you have built).

You can test your installation using the demonstration workspace
available at,

       http://www.network-theory.co.uk/gsl/demo.zip

Programs may need to be compiled with the following option selected,

       Project Settings
         C/C++
           Category: Customize
             Disable Language Extensions

depending on which header files you use.

The library is built with the /LD or /LDd option.  This is compatible
with the default link option /ML for single-threaded applications.
To use GSL in a multi-threaded application you may need to recompile
the library with another option, such as /MD.  See the Microsoft
Visual C++ Manual for details on link options.

The initial scripts were provided by José Miguel Buenaposada
(jmbuena@dia.fi.upm.es) and subsequently modified by Brian Gough
