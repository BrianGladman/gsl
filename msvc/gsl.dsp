# Microsoft Developer Studio Project File - Name="libgsl" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=libgsl - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "libgsl.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "libgsl.mak" CFG="libgsl - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "libgsl - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "libgsl - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "..\msvc" /I "." /I ".." /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ  /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od  /I "..\msvc" /I "." /I ".." /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ  /c
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 
# Begin Target

# Name "libgsl - Win32 Release"
# Name "libgsl - Win32 Debug"

# Begin Group "libgsl"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\version.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsl"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsl"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\complex_internal.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsl"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsl"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\templates_on.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsl"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsl"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\templates_off.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsl"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsl"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\gsl_math.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsl"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsl"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\gsl_pow_int.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsl"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsl"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\gsl_nan.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsl"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsl"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\gsl_machine.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsl"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsl"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\gsl_mode.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsl"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsl"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\gsl_precision.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsl"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsl"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslblas"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\blas\blas.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\blas\gsl_blas.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\blas\gsl_blas_types.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblas"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslblock"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\block\init.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\file.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\block.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\block_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\init_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\fprintf_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\fwrite_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_char.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_complex.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_complex_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_complex_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_complex_long_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_int.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_long.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_long_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_short.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_uchar.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_uint.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_ulong.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\block\gsl_block_ushort.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslblock"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslblock"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslcblas"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\cblas\sasum.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\saxpy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\scasum.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\scnrm2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\scopy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\sdot.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\sdsdot.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\sgbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\sgemm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\sgemv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\sger.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\snrm2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\srot.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\srotg.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\srotm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\srotmg.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ssbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\sscal.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\sspmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\sspr.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\sspr2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\sswap.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ssymm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ssymv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ssyr.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ssyr2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ssyr2k.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ssyrk.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\stbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\stbsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\stpmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\stpsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\strmm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\strmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\strsm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\strsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dasum.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\daxpy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dcopy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ddot.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dgbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dgemm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dgemv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dger.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dnrm2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\drot.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\drotg.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\drotm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\drotmg.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dsbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dscal.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dsdot.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dspmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dspr.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dspr2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dswap.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dsymm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dsymv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dsyr.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dsyr2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dsyr2k.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dsyrk.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dtbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dtbsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dtpmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dtpsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dtrmm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dtrmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dtrsm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dtrsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dzasum.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\dznrm2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\caxpy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ccopy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cdotc_sub.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cdotu_sub.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cgbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cgemm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cgemv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cgerc.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cgeru.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\chbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\chemm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\chemv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cher.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cher2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cher2k.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cherk.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\chpmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\chpr.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\chpr2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cscal.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\csscal.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cswap.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\csymm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\csyr2k.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\csyrk.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ctbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ctbsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ctpmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ctpsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ctrmm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ctrmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ctrsm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ctrsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zaxpy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zcopy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zdotc_sub.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zdotu_sub.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zdscal.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zgbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zgemm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zgemv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zgerc.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zgeru.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zhbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zhemm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zhemv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zher.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zher2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zher2k.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zherk.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zhpmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zhpr.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zhpr2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zscal.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zswap.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zsymm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zsyr2k.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\zsyrk.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ztbmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ztbsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ztpmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ztpsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ztrmm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ztrmv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ztrsm.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\ztrsv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\icamax.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\idamax.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\isamax.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\izamax.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\tests.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\tests.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\cblas.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_asum_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_asum_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_axpy_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_axpy_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_copy_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_copy_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_dot_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_dot_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_gbmv_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_gbmv_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_gemm_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_gemm_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_gemv_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_gemv_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_ger.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_gerc.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_geru.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_hbmv.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_hemm.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_hemv.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_her.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_her2.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_her2k.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_herk.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_hpmv.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_hpr.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_hpr2.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_iamax_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_iamax_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_nrm2_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_nrm2_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_rot.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_rotg.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_rotm.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_rotmg.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_sbmv.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_scal_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_scal_c_s.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_scal_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_spmv.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_spr.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_spr2.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_swap_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_swap_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_symm_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_symm_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_symv.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_syr.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_syr2.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_syr2k_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_syr2k_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_syrk_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_syrk_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_tbmv_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_tbmv_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_tbsv_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_tbsv_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_tpmv_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_tpmv_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_tpsv_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_tpsv_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_trmm_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_trmm_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_trmv_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_trmv_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_trsm_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_trsm_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_trsv_c.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\source_trsv_r.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\hypot.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cblas\gsl_cblas.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcblas"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcblas"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslcheb"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\cheb\deriv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcheb"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcheb"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cheb\eval.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcheb"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcheb"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cheb\init.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcheb"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcheb"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cheb\integ.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcheb"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcheb"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\cheb\gsl_chebyshev.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcheb"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcheb"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslcomplex"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\complex\math.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcomplex"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcomplex"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\complex\gsl_complex.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcomplex"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcomplex"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\complex\gsl_complex_math.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslcomplex"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslcomplex"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgsldht"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\dht\dht.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsldht"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsldht"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\dht\gsl_dht.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsldht"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsldht"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgsldiff"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\diff\diff.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsldiff"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsldiff"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\diff\gsl_diff.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsldiff"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsldiff"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgsleigen"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\eigen\eigen_sort.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsleigen"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsleigen"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\eigen\jacobi.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsleigen"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsleigen"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\eigen\gsl_eigen.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsleigen"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsleigen"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslerr"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\err\error.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslerr"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslerr"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\err\stream.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslerr"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslerr"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\err\message.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslerr"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslerr"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\err\strerror.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslerr"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslerr"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\err\warn.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslerr"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslerr"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\err\gsl_errno.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslerr"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslerr"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\err\gsl_message.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslerr"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslerr"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslfft"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\fft\dft.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\fft.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\c_pass.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\hc_pass.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\real_pass.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\signals.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\signals_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\c_main.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\c_init.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\c_pass_2.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\c_pass_3.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\c_pass_4.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\c_pass_5.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\c_pass_6.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\c_pass_7.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\c_pass_n.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\c_radix2.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\bitreverse.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\bitreverse.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\factorize.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\factorize.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\hc_init.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\hc_pass_2.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\hc_pass_3.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\hc_pass_4.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\hc_pass_5.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\hc_pass_n.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\hc_radix2.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\hc_unpack.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\real.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\real_init.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\real_pass_2.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\real_pass_3.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\real_pass_4.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\real_pass_5.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\real_pass_n.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\real_radix2.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\real_unpack.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\compare.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\compare_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\dft_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\hc_main.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\real_main.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\urand.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\gsl_fft.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\gsl_fft_complex.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\gsl_fft_halfcomplex.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\gsl_fft_real.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\gsl_dft_complex.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\gsl_dft_complex_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\gsl_fft_complex_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\gsl_fft_halfcomplex_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fft\gsl_fft_real_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfft"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfft"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslfit"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\fit\linear.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\fit\gsl_fit.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslfit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslfit"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslhistogram"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\histogram\add.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\get.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\init.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\params.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\reset.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\file.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\pdf.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\gsl_histogram.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\add2d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\get2d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\init2d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\params2d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\reset2d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\file2d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\pdf2d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\gsl_histogram2d.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\calloc_range.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\calloc_range2d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\copy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\copy2d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\maxval.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\maxval2d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\oper.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\oper2d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\find.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\find2d.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\gsl_histogram.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\histogram\gsl_histogram2d.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslhistogram"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslhistogram"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslieeeutils"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\ieee-utils\print.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\make_rep.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\gsl_ieee_utils.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\env.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\read.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-aix.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-darwin.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-hpux.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-hpux11.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-irix.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-m68klinux.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-ppclinux.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-solaris.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-sparclinux.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-sunos4.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-tru64.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-unknown.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-x86linux.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-freebsd.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-os2emx.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\fp-netbsd.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\endian.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\standardize.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ieee-utils\gsl_ieee_utils.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslieeeutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslieeeutils"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslintegration"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\integration\qk15.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qk21.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qk31.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qk41.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qk51.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qk61.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qk.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qng.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qng.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qag.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qags.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qagp.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\workspace.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qcheb.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qawc.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qmomo.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qaws.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qmomof.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qawo.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qawf.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qpsrt.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qpsrt2.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qelg.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qc25c.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qc25s.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\qc25f.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\ptsort.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\util.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\err.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\integration\gsl_integration.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslintegration"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslintegration"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslinterpolation"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\interpolation\accel.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslinterpolation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslinterpolation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\interpolation\akima.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslinterpolation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslinterpolation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\interpolation\bsearch.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslinterpolation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslinterpolation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\interpolation\cspline.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslinterpolation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslinterpolation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\interpolation\interp.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslinterpolation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslinterpolation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\interpolation\linear.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslinterpolation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslinterpolation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\interpolation\integ_eval_macro.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslinterpolation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslinterpolation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\interpolation\bsearch.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslinterpolation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslinterpolation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\interpolation\gsl_interp.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslinterpolation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslinterpolation"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgsllinalg"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\linalg\multiply.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\linalg\tridiag.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\linalg\tridiag.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\linalg\lu.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\linalg\hh.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\linalg\qr.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\linalg\qrpt.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\linalg\svd.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\linalg\householder.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\linalg\cholesky.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\linalg\givens.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\linalg\gsl_linalg.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsllinalg"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsllinalg"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslmatrix"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\matrix\init.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\matrix.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\file.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\rowcol.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\swap.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\copy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\minmax.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\prop.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\oper.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\getset.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\view.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\matrix_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\init_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\file_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\rowcol_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\swap_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\copy_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\minmax_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\prop_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\oper_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\getset_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\view_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_char.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_complex_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_complex_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_complex_long_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_int.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_long.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_long_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_short.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_uchar.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_uint.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_ulong.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\matrix\gsl_matrix_ushort.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmatrix"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmatrix"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslmin"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\min\fsolver.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\min\golden.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\min\brent.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\min\convergence.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\min\bracketing.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\min\min.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\min\gsl_min.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmin"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslmonte"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\monte\miser.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmonte"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmonte"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\monte\plain.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmonte"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmonte"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\monte\gsl_monte_vegas.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmonte"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmonte"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\monte\gsl_monte_miser.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmonte"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmonte"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\monte\gsl_monte_plain.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmonte"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmonte"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\monte\gsl_monte.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmonte"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmonte"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\monte\vegas.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmonte"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmonte"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\monte\gsl_monte.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmonte"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmonte"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\monte\gsl_monte_vegas.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmonte"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmonte"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\monte\gsl_monte_miser.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmonte"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmonte"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\monte\gsl_monte_plain.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmonte"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmonte"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslmultifit"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\multifit\multilinear.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\work.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\lmder.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\fsolver.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\fdfsolver.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\convergence.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\gradient.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\covar.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\lmutil.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\lmpar.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\lmset.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\lmiterate.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\qrsolv.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\gsl_multifit.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multifit\gsl_multifit_nlin.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultifit"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultifit"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslmultimin"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\multimin\directional_minimize.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultimin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultimin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multimin\fdfminimizer.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultimin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultimin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multimin\steepest_descent.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultimin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultimin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multimin\conjugate.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultimin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultimin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multimin\convergence.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultimin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultimin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multimin\vector_bfgs.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultimin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultimin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multimin\diff.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultimin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultimin"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multimin\gsl_multimin.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultimin"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultimin"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslmultiroots"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\multiroots\fdjac.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\fsolver.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\fdfsolver.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\convergence.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\newton.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\gnewton.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\dnewton.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\broyden.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\hybrid.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\hybridj.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\enorm.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\dogleg.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\multiroots\gsl_multiroots.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslmultiroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslmultiroots"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslntuple"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\ntuple\ntuple.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslntuple"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslntuple"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ntuple\gsl_ntuple.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslntuple"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslntuple"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslodeiv"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\ode-initval\control.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\evolve.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\monitor.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\bsimp.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\gear1.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\gear2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\odeiv.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\odeiv_util.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\rk2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\rk4.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\rk2imp.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\rk4imp.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\rk8pd.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\rkck.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\odeiv_util.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\ode-initval\gsl_odeiv.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslodeiv"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslodeiv"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslpermutation"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\permutation\init.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\file.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\permutation.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\permute.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\permute_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permutation.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_char.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_int.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_long.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_long_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_short.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_uchar.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_uint.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_ulong.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_ushort.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector_char.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector_int.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector_long.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector_long_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector_short.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector_uchar.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector_uint.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector_ulong.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\permutation\gsl_permute_vector_ushort.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpermutation"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpermutation"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslpoly"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\poly\eval.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\poly\solve_quadratic.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\poly\solve_cubic.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\poly\zsolve_quadratic.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\poly\zsolve_cubic.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\poly\zsolve.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\poly\zsolve_init.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\poly\balance.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\poly\companion.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\poly\norm.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\poly\qr.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\poly\gsl_poly.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslpoly"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslpoly"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslqrng"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\qrng\gsl_qrng.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslqrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslqrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\qrng\qrng.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslqrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslqrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\qrng\niederreiter-2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslqrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslqrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\qrng\sobol.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslqrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslqrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\qrng\gsl_qrng.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslqrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslqrng"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslrandist"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\randist\bernoulli.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\beta.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\bigauss.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\binomial.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\cauchy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\chisq.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\discrete.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\erlang.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\exponential.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\exppow.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\fdist.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\flat.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\gamma.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\gauss.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\gausstail.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\geometric.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\gumbel.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\hyperg.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\laplace.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\levy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\logarithmic.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\logistic.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\lognormal.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\nbinomial.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\pareto.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\pascal.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\poisson.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\rayleigh.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\shuffle.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\sphere.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\tdist.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\weibull.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\randist\gsl_randist.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrandist"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrandist"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslrng"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\rng\rng.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\types.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\default.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\cmrg.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\gfsr4.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\slatec.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\minstd.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\mrg.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\mt.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\r250.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\ran0.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\ran1.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\ran2.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\ran3.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\rand.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\random.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\rand48.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\randu.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\ranf.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\ranlux.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\ranlxs.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\ranlxd.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\ranmar.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\taus.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\transputer.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\tt.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\uni.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\uni32.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\vax.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\zuf.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\rng\gsl_rng.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslrng"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslrng"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslroots"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\roots\bisection.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\roots\brent.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\roots\falsepos.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\roots\newton.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\roots\secant.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\roots\steffenson.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\roots\convergence.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\roots\fsolver.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\roots\fdfsolver.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\roots\roots.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslroots"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\roots\gsl_roots.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslroots"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslroots"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslsiman"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\siman\siman.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsiman"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsiman"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\siman\gsl_siman.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsiman"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsiman"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslsort"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\sort\sort.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\sortind.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\sortvec.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\sortvecind.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\subset.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\subsetind.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\sortvec_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\sortvecind_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\subset_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\subsetind_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_heapsort.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_char.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_int.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_long.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_long_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_short.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_uchar.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_uint.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_ulong.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_ushort.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector_char.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector_int.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector_long.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector_long_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector_short.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector_uchar.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector_uint.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector_ulong.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sort\gsl_sort_vector_ushort.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsort"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsort"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslspecfunc"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\specfunc\airy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\airy_der.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\airy_zero.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\atanint.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_I0.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_I1.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_In.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_Inu.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_J0.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_J1.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_Jn.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_Jnu.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_K0.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_K1.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_Kn.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_Knu.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_Y0.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_Y1.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_Yn.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_Ynu.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_amp_phase.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_amp_phase.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_i.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_j.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_k.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_olver.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_temme.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_y.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_zero.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_sequence.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\beta.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\beta_inc.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\clausen.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\coulomb.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\coupling.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\coulomb_bound.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\dawson.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\debye.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\dilog.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\elementary.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\ellint.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\elljac.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\erfc.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\exp.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\expint.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\expint3.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\fermi_dirac.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gegenbauer.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gamma.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gamma_inc.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\hyperg_0F1.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\hyperg_2F0.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\hyperg_1F1.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\hyperg_2F1.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\hyperg_U.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\hyperg.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\laguerre.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\legendre_H3d.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\legendre_Qn.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\legendre_con.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\legendre_poly.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\log.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\poch.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\pow_int.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\psi.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\recurse.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\result.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\shint.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\sinint.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\synchrotron.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\transport.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\trig.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\zeta.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_amp_phase.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_olver.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel_temme.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\bessel.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\hyperg.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\legendre.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\eval.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\chebyshev.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\cheb_eval.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\cheb_eval_mode.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_specfunc.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_airy.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_bessel.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_clausen.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_coupling.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_coulomb.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_dawson.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_debye.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_dilog.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_elementary.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_ellint.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_elljac.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_erf.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_exp.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_expint.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_fermi_dirac.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_gegenbauer.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_gamma.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_hyperg.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_laguerre.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_legendre.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_log.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_pow_int.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_psi.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_result.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_synchrotron.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_transport.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_trig.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\specfunc\gsl_sf_zeta.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslspecfunc"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslspecfunc"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslstatistics"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\statistics\mean.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\variance.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\absdev.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\skew.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\kurtosis.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\lag1.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\p_variance.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\minmax.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\ttest.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\median.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\covariance.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\quantiles.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\wmean.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\wvariance.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\wabsdev.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\wskew.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\wkurtosis.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\mean_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\variance_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\covariance_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\absdev_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\skew_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\kurtosis_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\lag1_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\p_variance_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\minmax_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\ttest_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\median_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\quantiles_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\wmean_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\wvariance_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\wabsdev_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\wskew_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\wkurtosis_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics_char.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics_int.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics_long.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics_long_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics_short.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics_uchar.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics_uint.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics_ulong.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\statistics\gsl_statistics_ushort.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslstatistics"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslstatistics"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslsum"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\sum\levin_u.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsum"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsum"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sum\levin_utrunc.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsum"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsum"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sum\work_u.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsum"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsum"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sum\work_utrunc.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsum"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsum"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sum\gsl_sum.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsum"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsum"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslsys"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\sys\minmax.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsys"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsys"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sys\prec.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsys"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsys"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sys\hypot.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsys"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsys"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sys\log1p.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsys"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsys"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sys\expm1.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsys"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsys"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sys\coerce.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsys"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsys"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sys\invhyp.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsys"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsys"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\sys\pow_int.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslsys"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslsys"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgsltest"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\test\results.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsltest"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsltest"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\test\gsl_test.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgsltest"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgsltest"

!ENDIF 

# End Source File
# End Group
# Begin Group "libutils"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\utils\system.h



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libutils"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\utils\placeholder.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libutils"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libutils"

!ENDIF 

# End Source File
# End Group
# Begin Group "libgslvector"

# PROP Default_Filter ""

# Begin Source File

SOURCE=..\vector\init.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\file.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\vector.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\copy.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\swap.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\prop.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\minmax.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\oper.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\subvector.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\view.c



!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\vector_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\init_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\file_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\copy_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\swap_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\prop_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\minmax_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\oper_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\subvector_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\view_source.c

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_char.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_complex.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_complex_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_complex_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_complex_long_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_float.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_int.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_long.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_long_double.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_short.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_uchar.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_uint.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_ulong.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\vector\gsl_vector_ushort.h

# PROP Exclude_From_Build 1

!IF  "$(CFG)" == "libgsl - Win32 Release"

# PROP Intermediate_Dir "Release\libgslvector"

!ELSEIF  "$(CFG)" == "libgsl - Win32 Debug"

# PROP Intermediate_Dir "Debug\libgslvector"

!ENDIF 

# End Source File
# End Group
# End Target
# End Project
