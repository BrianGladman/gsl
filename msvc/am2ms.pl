#!/usr/bin/perl

# This script is a mess!

$target = shift @ARGV;

@confs = ("Release", "Debug");

if ($target eq 'GSLLIBML') {
    $begin_project_lib = \&begin_project_static_lib;
    @options = ("Release", "/ML", "Debug", "/MLd");
} elsif ($target eq 'GSLLIBMT') {
    $begin_project_lib = \&begin_project_static_lib;
    @options = ("Release", "/MT", "Debug", "/MTd");
#} elsif ($target eq 'GSLLIBMD') {
#    $begin_project_lib = \&begin_project_static_lib;
#    @options = ("Release", "/MD", "Debug", "/MDd");
} elsif ($target eq 'GSLDLL') {
    $begin_project_lib = \&begin_project_dll;
    @options = ("Release", "/MD /LD", "Debug", "/MDd /LDd");
    #$app_options = '/D "GSL_IMPORTS"';
} else {
    die "unrecognized target $target";
}


%automake = ();

for $file (@ARGV) {
    $/=undef;
    open(FILE,"<$file");
    $in = <FILE>;
    close(FILE);

    ($dir = $file) =~ s/\/Makefile.am$//;
    #$dir =~ s/.*\///;

    $in =~ s/\\\n/ /g;
    $in =~ s/\n\n/\n/g;
    $in =~ s/[ \t][ \t]+/ /g;
    $in =~ s/^\s+//g;
    
    {
        my %AM;
        @lines = split("\n",$in);
        for (@lines) {
            s/#.*//g;
            ($var,$data) = split(/\s*=\s*/, $_, 2);
            $AM{$var} = $data;
            # print "$var IS $data\n";
        }
        $automake{$file} = new Makefile ($dir, \%AM);
    }
}


for $file (@ARGV) {
    my $makefile = $automake{$file};

    @instlibs = $makefile->inst_libs();
    @noinstlibs = $makefile->noinst_libs();

    @gsllib = grep(/libgsl\.la/, @instlibs);
    @otherlibs = grep(!/libgsl\.la/, @instlibs);

    @ext_headers = $makefile->ext_headers();
    @int_headers = $makefile->int_headers();
    @progs = $makefile->programs();

    my $dir = $makefile->dir();

    for $lib (@otherlibs) {
        my $name = base($lib);
        push(@lib_projects, $name);
        #print "adding to $name \n";
        open(LIB, ">$name.dsp");
        print LIB &$begin_project_lib($name, @options);
        print LIB &begin_target($name);
        print LIB &add_def(".", "$name.def") if $target =~ /DLL/i;
        print LIB &add_files($name, @confs, $dir, $name, 'inc', $makefile->list_sources($lib));
        print LIB &add_files($name, @confs, $dir, $name, 'exclude', grep(!/^test_/,@int_headers));
        print LIB &add_files($name, @confs, $dir, $name, 'exclude', @ext_headers);
        print LIB &end_target();
        print LIB &end_project();
        close(LIB);
    }
}

open(LIBGSL, ">gsl.dsp");
print LIBGSL &$begin_project_lib("gsl", @options, "gslcblas.lib");
print LIBGSL &begin_target("gsl");
print LIBGSL &add_def(".", "gsl.def") if $target =~ /DLL/i;
push(@lib_projects, "gsl");
for $file (@ARGV) {
    my $makefile = $automake{$file};

    @instlibs = $makefile->inst_libs();
    @gsllib = grep(/libgsl\.la/, @instlibs);

    @noinstlibs = $makefile->noinst_libs();

    @ext_headers = $makefile->ext_headers();
    @int_headers = $makefile->int_headers();

    my $dir = $makefile->dir();

    push(@installed_headers, map("$dir/$_", @ext_headers));

    #print "== dir $dir ==\n" ;

    for $lib (@gsllib, @noinstlibs) {
        my $name = base($lib);
        #print "adding to libgsl ", $name, "\n";
        print LIBGSL &begin_group($name);
        print LIBGSL &add_files("gsl", @confs, $dir, $name, 'inc', $makefile->list_sources($lib));
        print LIBGSL &add_files("gsl", @confs, $dir, $name, 'exclude', grep(!/^test_/,@int_headers));
        print LIBGSL &add_files("gsl", @confs, $dir, $name, 'exclude', @ext_headers);
        print LIBGSL &end_group();
    }
}
print LIBGSL &end_target();
print LIBGSL &end_project();

for $file (@ARGV) {
    my $makefile = $automake{$file};

    @instlibs = $makefile->inst_libs();
    @noinstlibs = $makefile->noinst_libs();

    @gsllib = grep(/libgsl\.la/, @instlibs);
    @otherlibs = grep(!/libgsl\.la/, @instlibs);

    @ext_headers = $makefile->ext_headers();
    @int_headers = $makefile->int_headers();
    @progs = $makefile->programs();

    my $dir = $makefile->dir();

    for $prog (@progs) {
        #print "put in tests $prog\n";
        @sources = $makefile->list_sources($prog);
        @ldadds = $makefile->list_ldadds($prog);
        my $name = $dir; 
        $name =~ s/.*\///; 
        $name =~ s/-/_/g;
        if ($prog =~ /^test/) {
            $name = "${prog}_${name}";
            push (@tests, "$name");
        } else {
            $name = $prog;
        }
        open(TEST, ">$name.dsp");
        warn "$name $prog => ", @sources, "\n";
        print TEST &begin_project_app($name, @options, $app_options);
        print TEST &begin_target($name);
        print TEST &add_files($name, @confs, $dir, $name, 'inc', $makefile->list_sources($prog));
        print TEST &add_files($name, @confs, $dir, $name, 'exclude', grep(/^test_/,@int_headers));
        print TEST &end_target();
        print TEST &end_project();
        close (TEST);
    }
}

#  Write GSL workspace

open(GSL, ">GSL.dsw");
print GSL &begin_workspace();
for $lib (@lib_projects) {
    print GSL &add_workspace_project("$lib", "$target\\$lib.dsp");
}
print GSL &end_workspace();
close GSL;

open(GSL, ">GSLTESTS.dsw");
print GSL &begin_workspace();
for $t (@tests) {
    print GSL &add_workspace_project($t, "$target\\$t.dsp");
}
print GSL &end_workspace();
close GSL;

# Write test batch files

for $dir (@confs) {
    open(TEST, ">MAKE_CHECK_${dir}.bat");
    print TEST &begin_check();

    for $t (@tests) {
        print TEST &add_check($dir,$t);
    }

    print TEST &missing_checks();
    print TEST &end_check();
    close(TEST);
}

######################################################################

sub begin_project_static_lib {
    my ($proj, $release, $release_options, $debug, $debug_options) = @_;

    return <<"EOF";
# Microsoft Developer Studio Project File - Name="$proj" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=$proj - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "$proj.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "$proj.mak" CFG="$proj - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "$proj - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "$proj - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "\$(CFG)" == "$proj - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "$release"
# PROP BASE Intermediate_Dir "$release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "$release"
# PROP Intermediate_Dir "$release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /Op- /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo ${release_options} /Za /W3 /GX /O2 /Op- /I "..\\msvc" /I "..\\gsl" /I "." /I ".." /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "\$(CFG)" == "$proj - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "$debug"
# PROP BASE Intermediate_Dir "$debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "$debug"
# PROP Intermediate_Dir "$debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ  /c
# ADD CPP /nologo ${debug_options} /Za /W3 /Gm /GX /Z7 /Od  /I "..\\msvc" /I "..\\gsl" /I "." /I ".." /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ  /c
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 
EOF
}

sub begin_project_dll {
    my ($proj, $release, $release_options, $debug, $debug_options, $extra_libs) = @_;

    return <<"EOF";
# Microsoft Developer Studio Project File - Name="$proj" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=$proj - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "$proj.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "$proj.mak" CFG="$proj - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "$proj - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "$proj - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "\$(CFG)" == "$proj - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "$release"
# PROP BASE Intermediate_Dir "$release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "$release"
# PROP Intermediate_Dir "$release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo ${release_options} /W3 /GX /O2 /Op- /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "GSL_DLL" /D "DLL_EXPORT" /YX /FD /c
# ADD CPP /nologo  ${release_options} /Za /W3 /GX /O2 /Op- /I "..\\msvc" /I "..\\gsl" /I "." /I ".." /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "GSL_DLL" /D "DLL_EXPORT" /YX /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386
# ADD LINK32 $extra_libs kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386 /libpath:"$release"

!ELSEIF  "\$(CFG)" == "$proj - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "$debug"
# PROP BASE Intermediate_Dir "$debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "$debug"
# PROP Intermediate_Dir "$debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo  ${debug_options} /W3 /Gm /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "GSL_DLL" /D "DLL_EXPORT" /YX /FD /GZ /c
# ADD CPP /nologo  ${debug_options} /Za /W3 /Gm /GX /Z7 /Od /I "..\\msvc" /I "..\\gsl" /I "." /I ".." /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "GSL_DLL" /D "DLL_EXPORT" /FD /LDd /GZ /c
# SUBTRACT CPP /YX
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /fo"$proj.res" /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 $extra_libs kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /verbose /dll /incremental:no /debug /machine:I386 /pdbtype:sept /libpath:"$debug"

!ENDIF 
EOF
}


sub begin_project_app {
    my ($proj, $release, $release_options, $debug, $debug_options, $app_options) = @_;

    return <<"EOF";
# Microsoft Developer Studio Project File - Name="$proj" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=$proj - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "$proj.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "$proj.mak" CFG="$proj - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "$proj - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "$proj - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "\$(CFG)" == "$proj - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "$release"
# PROP BASE Intermediate_Dir "$release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "$release"
# PROP Intermediate_Dir "$release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /Op- /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo ${release_options} /Za /W3 /GX /O2 /Op- /I "..\\msvc" /I "..\\gsl" /I "." /I ".." /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" ${app_options} /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib  kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 gsl.lib gslcblas.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386 /out:"$release\\$proj.exe" /libpath:"$release"

!ELSEIF  "\$(CFG)" == "$proj - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "$debug"
# PROP BASE Intermediate_Dir "$debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "$debug"
# PROP Intermediate_Dir "$debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ  /c
# ADD CPP /nologo ${debug_options} /Za /W3 /Gm /GX /Z7 /Od /I "..\\msvc" /I "..\\gsl" /I "." /I ".." /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" ${app_options} /FD /GZ  /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib  kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 gsl.lib gslcblas.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /out:"$debug\\$proj.exe" /pdbtype:sept /libpath:"$debug"

!ENDIF 
EOF
}


sub end_project {
    return <<'EOF';
# End Project
EOF
}

#----------------------------------------------------------------------

sub begin_target {
    my ($proj) = @_;
    return <<"EOF";

# Begin Target

# Name "$proj - Win32 Release"
# Name "$proj - Win32 Debug"
EOF
}

sub end_target {
    return <<'EOF';
# End Target
EOF
}

#----------------------------------------------------------------------

sub begin_group {
    my ($group) = @_;
    return <<"EOF";
# Begin Group "$group"

# PROP Default_Filter ""

EOF
}

sub end_group {
    return <<'EOF';
# End Group
EOF
}

#----------------------------------------------------------------------

sub add_def {
    my ($dir, $file) = @_;
    my $y;
    $dir =~ s/\//\\/g; # convert to dos style path with \'s
    return <<"EOF";
# Begin Source File

SOURCE=$dir\\$file
# End Source File
EOF
}

sub add_files {
    my ($proj, $release, $debug, $dir, $name, $flag, @files) = @_;
    my $file;
    my $y;
    $dir =~ s/\//\\/g; # convert to dos style path with \'s
    my $prop;
    $prop = "# PROP Exclude_From_Build 1" if $flag eq 'exclude';
    for $file (@files) {
        my $location;
        $file =~ s/\//\\/g; # convert to dos style path with \'s
        if ($file =~ /gsl_/) {
            $location = $dir;
            $location =~ s#\\[^\\]*$##;
            $location .= "\\gsl\\$file";
        } else {
            $location = "$dir\\$file";
        }
    my $x = 
"# Begin Source File

SOURCE=$location

$prop

!IF  \"\$(CFG)\" == \"$proj - Win32 Release\"

# PROP Intermediate_Dir \"$release\\$name\"

!ELSEIF  \"\$(CFG)\" == \"$proj - Win32 Debug\"

# PROP Intermediate_Dir \"$debug\\$name\"

!ENDIF 

# End Source File
";
    $y .= $x;}
    return $y;
}

######################################################################

sub begin_workspace {
    return <<'EOF';
Microsoft Developer Studio Workspace File, Format Version 6.00
# WARNING: DO NOT EDIT OR DELETE THIS WORKSPACE FILE!

EOF
}

sub end_workspace {
    return <<'EOF';
###############################################################################

Global:

Package=<5>
{{{
}}}

Package=<3>
{{{
}}}

###############################################################################

EOF
}

sub add_workspace_project {
    my ($name, $file, @dependencies) = @_;

    my $depends = "";
    for my $d (@dependencies) {
        $depends .= "    Begin Project Dependency\n";
        $depends .= "    Project_Dep_Name $d\n";
        $depends .= "    End Project Dependency\n";
    }

    return <<"EOF";
###############################################################################

Project: "$name"="$file" - Package Owner=<4>

Package=<5>
{{{
}}}

Package=<4>
{{{
$depends}}}

EOF
}

######################################################################

sub generic_project {
    my ($name, $release, $debug) = @_;
    return <<"EOF";
# Microsoft Developer Studio Project File - Name="$name" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Generic Project" 0x010a

CFG=$name - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "$name.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "$name.mak" CFG="$name - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "$name - Win32 Release" (based on "Win32 (x86) Generic Project")
!MESSAGE "$name - Win32 Debug" (based on "Win32 (x86) Generic Project")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
MTL=midl.exe

!IF  "\$(CFG)" == "$name - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "$release"
# PROP BASE Intermediate_Dir "$release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "$release"
# PROP Intermediate_Dir "$release"
# PROP Target_Dir ""

!ELSEIF  "\$(CFG)" == "$name - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "$debug"
# PROP BASE Intermediate_Dir "$debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "$debug"
# PROP Intermediate_Dir "$debug"
# PROP Target_Dir ""

!ENDIF 

# Begin Target

# Name "$name - Win32 Release"
# Name "$name - Win32 Debug"
# End Target
# End Project
EOF
}

######################################################################

sub begin_check {
    return <<'EOF';
@echo off
@REM
@REM  MS-DOS Batch file to test the GNU Scientific Library (GSL)
@REM  on M$ Visual C++ 6.0.
@REM
@REM  2001/04/03 - José Miguel Buenaposada (jmbuena@dia.fi.upm.es)
@REM               Initial version.
@REM  
@REM

@echo "Testing GSL library" 
if exist result.dat del result.dat
EOF
} 

sub add_check {
    my ($dir, $test) = @_;
    return <<"EOF";
\@echo running $dir $test...
\@$dir\\$test.exe >> result.dat
\@if not errorlevel 1 echo PASS: $dir\\$test >> result.dat
\@if errorlevel 1 echo FAIL: $dir\\$test >> result.dat
EOF
}

sub missing_checks {
    return <<'EOF';
@echo PASS: test_gsl_histogram.sh NOT IMPLEMENTED ON WINDOWS >> result.dat
EOF
}

sub end_check {
    return <<'EOF'
@echo Number of failures:
@FIND /V /C "PASS:" result.dat
@echo Number of passes:
@FIND /C "PASS:" result.dat
@echo Full test log is in result.dat
EOF
}

######################################################################


sub base {
    my ($f) = @_ ;
    $f =~ s/\.\w+$//;
    $f =~ s/^lib//;
    return $f;
}

package Makefile;

sub new {
    my ($class, $dir, $args) = @_;
    my $self = {} ;
    bless $self;
    $self->{'dir'} = $dir;
    my %am = %$args;
    for my $key (keys %am) {
        $self->{AM}->{$key} = [split(' ', $am{$key})];
    }
    return $self;
}

sub dir {
    my ($self) = @_;
    return $self->{dir};
}


sub list_sources {
    my ($self, $f) = @_;
    $f =~ s/[\.\-]/_/g;
    return @{$self->{'AM'}->{"${f}_SOURCES"}};
}

sub list_ldadds {
    my ($self, $f) = @_;
    $f =~ s/[\.\-]/_/g;
    return @{$self->{'AM'}->{"${f}_LDADD"}};
}

sub inst_libs {
    my ($self) = @_;
    return @{$self->{'AM'}->{lib_LTLIBRARIES}};
}

sub noinst_libs {
    my ($self) = @_;
    return @{$self->{'AM'}->{noinst_LTLIBRARIES}};
}

sub ext_headers {
    my ($self) = @_;
    return @{$self->{'AM'}->{pkginclude_HEADERS}};
}

sub int_headers {
    my ($self) = @_;
    return @{$self->{'AM'}->{noinst_HEADERS}};
}

sub programs {
    my ($self) = @_;
    return (@{$self->{'AM'}{bin_PROGRAMS}}, @{$self->{'AM'}->{check_PROGRAMS}});
}

######################################################################
