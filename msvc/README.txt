This file contains instructions for compiling applications with GSL --
you may want to print it out for reference.

The GNU Scientific Library (for Microsoft Visual C++)
=====================================================

This is GSL, the GNU Scientific Library, a collection of numerical
routines for scientific computing.

GSL is free software, you can redistribute it and/or modify it under
the terms of the GNU General Public License.

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

More information about GSL
==========================

The project homepage is http://www.gnu.org/software/gsl/

The developers page is http://sources.redhat.com/gsl/

Please report bugs to the GSL discussion list gsl-discuss@sources.redhat.com

See the NEWS file for recent developments.

The GSL Manual has been published and can be ordered from most
bookstores. The publication details are,

  GNU Scientific Library Reference Manual, M. Galassi et al, ISBN
  095416170X (600 pages, paperback).

If you are interested in participating in GSL development, please send
mail to Mark Galassi -- rosalia@lanl.gov and Brian Gough --
bjg@network-theory.co.uk.

Using GSL with Microsoft Visual C++ 6.0
=======================================

The following sections describe how to compile an application with
GSL using Microsoft Visual C++.

By default the GSL libraries and header files are installed in the
following locations,

  C:\Program Files\GSL\include\gsl          - header files

  C:\Program Files\GSL\lib\libgsl.lib       - Main lib file
  C:\Program Files\GSL\lib\libgslcblas.lib  - BLAS lib file

  C:\Windows\System\libgsl.dll              - Main DLL
  C:\Windows\System\libgslcblas.dll         - BLAS DLL

Only the DLL 'Release' version of the library is supplied with this
package.  If you need other versions of the library you can build them
by recompiling the *.dsw project workspace files in the top-level src/
directory.

Compiling an Application
========================

To compile an application you will need to specify locations of the
GSL include files and libraries.  The installation program will add
the following global settings to the Windows Registry automatically,

    Tools --  Options -- Directories
          Include Files
             "C:\Program Files\GSL\include"
          Library Files: 
             "C:\Program Files\GSL\lib"

Alternatively they can be set for each project,

    Project Settings -- C/C++
        Category: Preprocessor  
          Additional Include Directories: 
             "C:\Program Files\GSL\include"

    Project Settings -- Link
        Category: Input
          Additional Library Path: 
             "C:\Program Files\GSL\lib"

You will always need to add the GSL .lib files to the list of
libraries for each project,

    Project Settings -- Link
        Category: Input
          Object/Library Modules: 
             libgsl.lib libgslcblas.lib 

Make sure that the Object/Library module settings are made for all the
appropriate configurations (either 'Release' or 'Debug').

The supplied DLLs are compiled with /MD, which generates thread-safe
Release DLLs.  You will need to compile your code with /MD in order to
link to the DLL.  The /MD option also defines the preprocessor
variable _DLL, which ensures that the appropriate DLL functions are
imported in the header files.

If you want to use the inline functions from GSL, you should also add
the preprocessor definitions HAVE_INLINE,inline=__inline.  The inline
functions are faster, but increase the code size.

You can test your installation using the demonstration workspace
available in the directory 'demo'.

Your programs will also need to be compiled with the following option
selected,

    Project Settings
      C/C++
        Category: Customize
          Disable Language Extensions

Sometimes this is not necessary, it depends on which header files your
application includes.  If you do apply this option you should use the
optimization option /Op- to prevent your code from being slowed down
by 'Disable Language extensions'.  By default the 'Disable Language
Extensions' option turns on the strictest IEEE arithmetic behavior,
which slows down the program significantly.  For most programs this is
not required and can be turned off with /Op-.

Single-Threaded vs Multi-Threaded Libraries
===========================================

Single-threaded libraries can be built from source with the /ML
option.  This is compatible with the default link option /ML for
single-threaded applications in Microsoft Visual C++.  Statically
linked multi-threaded libraries can be built with the /MT option.
Different versions of the library using these options can be built
from the workspaces GSLLIBML.dsw and GSLLIBMT.dsw.  See the Microsoft
Visual C++ Manual for more details on link options.

Building GSL with different compilation options
===============================================

There are three project workspaces for compiling the library as a DLL
or with the /ML (static single-threaded) and /MT (static
multi-threaded) options.  To compile the library use the following
procedure:

1.  Build the workspace GSLDLL.dsw, GSLLIBML.dsw, or GSLLIBMT.dsw
using 'Batch Build'.  These should generate their output in the
GSLDLL, GSLLIBML and GSLLIBMT directories, producing the following
libraries:

    GSLDLL/Debug/libgsl.{exp,lib,dll}
    GSLDLL/Debug/libgslcblas.{exp,lib,dll}
    GSLDLL/Release/libgsl.{exp,lib,dll}
    GSLDLL/Release/libgslcblas.{exp,lib,dll}

    GSLLIBML/Debug/libgsl.lib
    GSLLIBML/Debug/libgslcblas.lib
    GSLLIBML/Release/libgsl.lib
    GSLLIBML/Release/libgslcblas.lib

    GSLLIBMT/Debug/libgsl.lib
    GSLLIBMT/Debug/libgslcblas.lib
    GSLLIBMT/Release/libgsl.lib
    GSLLIBMT/Release/libgslcblas.lib

The header file msvc/config.h contains the appropriate configuration
for Microsoft Visual C++ and overrides the normal config.h.

2.  Build the test workspace GSLDLLTESTS.dsw, GSLLIBMLTESTS.dsw,
GSLLIBMTTESTS.dsw

This should create the test programs in GSLDLL, GSLLIBML, GSLLIBMT
directories.  To run the tests use the batch files,

        MAKE_CHECK_Debug.bat
        MAKE_CHECK_Release.bat

in each of these directories.  The MAKE CHECK scripts produce an
output log file "results.dat".  Any lines which don't begin with PASS:
indicate a problem.

3.  If the tests are successful you can link against the files
libgsl.lib, libgslcblas.lib and gsl/gsl*.h.  The *.dll files are also
needed to use the GSLDLL build.  You will need to add the path to
these files to the project settings for programs that you compile with
GSL.

Acknowledgements
================

The initial scripts to generate the Visual Studio project files for
GSL were provided by José Miguel Buenaposada (jmbuena@dia.fi.upm.es)
and subsequently modified by Brian Gough

GSL has been written by the GSL Team.  See the file AUTHORS.txt for
details.

