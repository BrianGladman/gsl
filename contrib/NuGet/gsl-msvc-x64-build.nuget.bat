REM @echo off

Echo LIB Windows Build NuGet

REM # XEON x64 Build Vars #
set _SCRIPT_DRIVE=%~d0
set _SCRIPT_FOLDER=%~dp0
set INITDIR=%_SCRIPT_FOLDER%
set SRC=%_SCRIPT_FOLDER%\..\..\
set BUILDTREE=%SRC%\build-win\
SET tbs_arch=x64
SET vcvar_arg=x86_amd64

REM # VC Vars #
SET VCVAR="%programfiles(x86)%\Microsoft Visual Studio\2017\Professional\VC\Auxiliary\Build\vcvarsall.bat"
if exist %VCVAR% call %VCVAR% %vcvar_arg%
SET VCVAR="%programfiles(x86)%\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvarsall.bat"
if exist %VCVAR% call %VCVAR% %vcvar_arg%

REM # Clean Build Tree #
rd /s /q %BUILDTREE%
mkdir %BUILDTREE%
cd %BUILDTREE%

SET platform=%tbs_arch%
set BINDIR=%SRC%\build-nuget\
rd /s /q %BINDIR%
mkdir %BINDIR%
echo %BINDIR%

:Release
msbuild /m %SRC%build.vc\gsl.lib.sln /p:Configuration=Release /p:Platform=%platform%
msbuild /m %SRC%build.vc\gsl.dll.sln /p:Configuration=Release /p:Platform=%platform%

:Debug
msbuild /m %SRC%build.vc\gsl.lib.sln /p:Configuration=Debug /p:Platform=%platform%
msbuild /m %SRC%build.vc\gsl.dll.sln /p:Configuration=Debug /p:Platform=%platform%

:copy
mkdir %BINDIR%lib\%platform%
mkdir %BINDIR%dll\%platform%
xcopy /S /I %SRC%lib\%platform% %BINDIR%lib\%platform%
xcopy /S /I %SRC%dll\%platform% %BINDIR%dll\%platform%

del %BINDIR%gsl
xcopy /I %SRC%gsl %BINDIR%gsl
copy %INITDIR%\gsl-msvc-%tbs_arch%.targets %BINDIR%\gsl-msvc-%tbs_arch%.targets

:nuget_req
REM # make nuget packages from binaries #
cd %BUILDTREE%
nuget pack %INITDIR%\gsl-msvc-%tbs_arch%.nuspec
cd %INITDIR%
REM --- exit ----
GOTO:eof
