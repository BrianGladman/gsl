This file contains instructions for compiling applications with GSL --
you may want to print it out for reference.

The GNU Scientific Library (Windows Version)
============================================

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

The project homepage is http://sources.redhat.com/gsl/

Please report bugs to the GSL discussion list gsl-discuss@sources.redhat.com

See the NEWS file for recent developments.

The GSL Manual has been published and can be ordered from most
bookstores. The publication details are,

  GNU Scientific Library Reference Manual, M. Galassi et al, ISBN
  095416170X (600 pages, paperback).

If you are interested in participating in GSL development, please send
mail to Mark Galassi -- rosalia@lanl.gov

Using GSL with Microsoft Visual C++ 6.0
=======================================

The following sections describe how to compile an application with
GSL using Microsoft Visual C++.

By default the GSL libraries and header files are installed in the
following locations,

    C:\Program Files\GSL\include\gsl  - header files
    C:\Program Files\GSL\lib          - lib files
    C:\Windows\System                 - DLL files

The "Release" versions of the libraries are installed in the lib
directory as follows,

    libgsl.lib,    libgslcblas.lib    - Release (DLL)
    libgslML.lib,  libgslcblasML.lib  - Release (single-threaded /ML)
    libgslMT.lib,  libgslcblasMT.lib  - Release (multi-threaded /MT)

If you need Debug versions of the library you can build them by
compiling the *.dsw project workspace files in the top-level src/
directory.

Compiling an Application
========================

To compile an application you will need to specify locations of the
GSL include files and libraries.  The installation program will add
the following global settings automatically,

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
             libgsl.lib libgslcblas.lib ... for the release DLL configuration
          or libgslMT.lib libgslcblasMT.lib ... (as above, multi-threaded)

Make sure that the Object/Library module settings are made for all the
appropriate configurations (either 'Release' or 'Debug').

If you are using the DLL you will need to set a preprocessor
definition for GSL_IMPORTS to import the DLL functions correctly.

If you want to use the inline functions from GSL, you should also add
the proprocessor definitions HAVE_INLINE,inline=__inline.  The inline
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

The single-threaded libraries are built with the /ML option.  This is
compatible with the default link option /ML for single-threaded
applications in Microsoft Visual C++.  To use GSL in a multi-threaded
application you will need to use the multi-threaded versions of the
library and compile your code with either /MT for the statically
linked version or /MD for the DLL version.  See the Microsoft Visual
C++ Manual for details on link options.

Acknowledgements
================

The initial scripts to generate the Visual Studio project files for
GSL were provided by José Miguel Buenaposada (jmbuena@dia.fi.upm.es)
and subsequently modified by Brian Gough
