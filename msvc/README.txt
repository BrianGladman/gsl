This file answers all the most common questions about
compiling applications with GSL -- you may want to
print it out for reference.

If the use of compiler options described here is
unfamiliar to you, refer to the online help supplied
with the Microsoft compiler itself for more
information.

The GNU Scientific Library (for Microsoft Visual C++)
=====================================================

This is GSL, the GNU Scientific Library, a collection
of numerical routines for scientific computing.

GSL is free software, you can redistribute it and/or
modify it under the terms of the GNU General Public
License.

The GNU General Public License does not permit this
software to be redistributed in proprietary programs.

This library is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

What is GNU?
============

GSL is part of the GNU Project.  The GNU Project was
launched in 1984 to develop a complete Unix-like
operating system which is free software: the GNU
system.  Variants of the GNU operating system are now
widely used.

The GNU Project is not limited to operating systems.
It aims to provide a whole spectrum of free software,
whatever many users want to have.

What is Free Software?
======================

Free Software is about freedom. It refers to the
freedom to share and modify software. The important
thing about free software is that everyone has the
freedom to cooperate with others in using and
improving it.

More information about GSL
==========================

The project homepage is http://www.gnu.org/software/gsl/

See the NEWS file for recent developments.

The GSL Manual has been published and can be ordered
from most bookstores. The publication details are,

  GNU Scientific Library Reference Manual, M. Galassi
  et al, ISBN 095416170X (600 pages, paperback).

If you are interested in participating in GSL
development, please send mail to
gsl-discuss@sources.redhat.com.

Using GSL with Microsoft Visual C++ 6.0
=======================================

The following sections describe how to compile an
application with GSL using Microsoft Visual C++.

By default the GSL libraries and header files are
installed in the following locations,

  C:\Program Files\GSL\include\gsl          - header files

  C:\Program Files\GSL\lib\libgsl.lib       - Main lib file
  C:\Program Files\GSL\lib\libgslcblas.lib  - BLAS lib file

Only the single-threaded 'Release' version of the
library is supplied with this package.  If you need
other versions of the library you can build them by
recompiling the *.dsw project workspace files in the
top-level src/ directory.

Compiling an Application
========================

To compile an application you will need to specify
locations of the GSL include files and libraries.
The installation program will add the following
global settings to the Windows Registry
automatically,

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

You will always need to add the GSL .lib files to the
list of libraries for each project,

    Project Settings -- Link
        Category: Input
          Object/Library Modules: 
             libgsl.lib libgslcblas.lib 

Make sure that the Object/Library module settings are
made for all the appropriate configurations (either
'Release' or 'Debug').

The supplied libraries are compiled with /ML, which
generates single-threaded release libraries.

If you want to use the inline functions from GSL, you
should also add the following preprocessor
definitions:

    HAVE_INLINE,inline=__inline

The inline functions are faster, but increase the
code size.

You can test your installation using the
demonstration workspace available in the directory
'demo'.

Your programs will also need to be compiled with the
following option selected,

    Project Settings
      C/C++
        Category: Customize
          Disable Language Extensions

Sometimes this is not necessary, it depends on which
header files your application includes.  If you do
apply this option you should use the optimization
option /Op- to prevent your code from being slowed
down by 'Disable Language extensions'.  By default
the 'Disable Language Extensions' option turns on the
strictest IEEE arithmetic behavior, which slows down
the program significantly.  For most programs this is
not required and can be turned off with /Op-.

Single-Threaded vs Multi-Threaded Libraries
===========================================

Statically linked single-threaded libraries are
supplied by default and are built from the source
with the /ML option.  This is compatible with the
default link option /ML for single-threaded
applications in Microsoft Visual C++.

Statically linked multi-threaded libraries can be
built with the /MT option.  Different versions of the
library using this options can be built from the
workspace GSLLIBMT.dsw.

Dynamically linked multi-threaded libraries can be
built with the /MD option.  Different versions of the
library using this options can be built from the
workspace GSLLIBMT.dsw.  You will need to compile
your code with /MD in order to link to the DLL.  The
/MD option can be enabled using the following
setting:

    Project Settings
      C/C++
        Category: Code Generation
          Multithreaded DLL

The /MD option also defines the preprocessor variable
_DLL, which ensures that the appropriate DLL
functions are imported in the header files.

To run the program compiled with 'ReleaseDLL' option
you will need to make sure the DLL files libgsl.dll
and libgslcblas.dll are available in the
C:\Windows\System directory, or copy them into the
same directory as the executable.

See the Microsoft Visual C++ Manual for more details
on link options.

Building GSL with different compilation options
===============================================

There are three project workspaces for compiling the
library with the /ML option (static single-threaded),
with the /MT option (static multi-threaded), or as a
DLL.  To compile the library use the following
procedure:

1.  Build the workspace GSLLIBML.dsw, GSLLIBMT.dsw or
GSLLIBMT.dsw using 'Batch Build'.  These should
generate their output in the GSLLIBMT, GSLLIBML and
GSLDLL directories, producing the following
libraries:

    GSLLIBML/Debug/libgsl.lib
    GSLLIBML/Debug/libgslcblas.lib
    GSLLIBML/Release/libgsl.lib
    GSLLIBML/Release/libgslcblas.lib

    GSLLIBMT/Debug/libgsl.lib
    GSLLIBMT/Debug/libgslcblas.lib
    GSLLIBMT/Release/libgsl.lib
    GSLLIBMT/Release/libgslcblas.lib

    GSLDLL/Debug/libgsl.{exp,lib,dll}
    GSLDLL/Debug/libgslcblas.{exp,lib,dll}
    GSLDLL/Release/libgsl.{exp,lib,dll}
    GSLDLL/Release/libgslcblas.{exp,lib,dll}

The header file msvc/config.h contains the
appropriate configuration for Microsoft Visual C++
and overrides the normal config.h.

2.  Build the test workspaces GSLLIBMLTESTS.dsw,
GSLLIBMTTESTS.dsw and GSLDLLTESTS.dsw.

This should create the test programs in GSLLIBML,
GSLLIBMT, GSLDLL directories.  To run the tests use
the batch files,

        MAKE_CHECK_Debug.bat
        MAKE_CHECK_Release.bat

in each of these directories.  The MAKE CHECK scripts
produce an output log file "results.dat".  Any lines
which don't begin with PASS: indicate a problem.

3.  If the tests are successful you can link against
the files libgsl.lib, libgslcblas.lib and gsl/gsl*.h.
The *.dll files are also needed to use the GSLDLL
build.  You will need to add the path to these files
to the project settings for programs that you compile
with GSL.

Reporting Bugs
==============

A list of known bugs can be found in the BUGS file.
Details of compilation problems can be found in the
INSTALL file.

If you find a bug which is not listed in these files
please report it to bug-gsl@gnu.org.

All bug reports should include:

 The version number of GSL, and where you obtained it.
 The hardware and operating system
 The compiler used, including version number and compilation options
 A description of the bug behaviour
 A short program which reproducibly exercises the bug

Thank you.

Acknowledgements
================

The initial scripts to generate the Visual Studio
project files for GSL were provided by José Miguel
Buenaposada (jmbuena@dia.fi.upm.es) and subsequently
modified by Brian Gough

GSL has been written by the GSL Team.  See the file
AUTHORS.txt for details.

