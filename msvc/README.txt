The GNU Scientific Library (Windows Version)
============================================

This is GSL, the GNU Scientific Library.

GSL is free software, you can redistribute it and/or modify it under
the terms of the GNU General Public License.

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

This file contains instructions for compiling applications with GSL --
you may want to print it out for reference.

More information about GSL
==========================

The library is in the final stages of development, and is being made
available in hopes that others might be interested in joining the
project and debugging the code.

See the NEWS file for recent developments.

The project homepage is http://sources.redhat.com/gsl/

Please report bugs to the GSL discussion list gsl-discuss@sources.redhat.com

If you are interested in participating in GSL development, please send
mail to Mark Galassi -- rosalia@lanl.gov

Using GSL with Microsoft Visual C++ 6.0
=======================================

By default the GSL libraries and header files are installed in the
following locations,

    C:\Program Files\GSL\include\gsl  - header files
    C:\Program Files\GSL\lib          - lib files

The "Release" and "Debug" versions of the libraries are installed in
the lib directory as follows,

    libgsl.lib,  libgslcblas.lib   - release version
    libgsld.lib, libgslcblasd.lib  - debug version

The debug version has a 'd' at the end of the library name.

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
             libgsl.lib libgslcblas.lib ...   for the release configuration
          or libgsld.lib libgslcblasd.lib ... for the debug configuration

Make sure that the Object/Library module settings are made for all the
appropriate configurations (either 'Release' or 'Debug').

If you want to use the inline functions from GSL, you should also add
the proprocessor definitions HAVE_INLINE,inline=__inline.  The inline
functions are faster, but increase the code size.

You can test your installation using the demonstration workspace
available in the directory 'demo'.

Depending on which header files your application uses some programs
may need to be compiled with the following option selected,

    Project Settings
      C/C++
        Category: Customize
          Disable Language Extensions

If you do this you should apply the optimization option /Op- to
prevent your code from being slowed down by a feature of 'Disable
Language extensions'.  By default the 'Disable Language Extensions'
option turns on the strictest IEEE arithmetic behavior, which slows
down the program significantly.  For most programs this is not
required and can be turned off with /Op-.

The library is built with the /LD or /LDd option.  This is compatible
with the default link option /ML for single-threaded applications.
To use GSL in a multi-threaded application you may need to recompile
the library with another option, such as /MD.  See the Microsoft
Visual C++ Manual for details on link options.

The initial scripts were provided by José Miguel Buenaposada
(jmbuena@dia.fi.upm.es) and subsequently modified by Brian Gough
