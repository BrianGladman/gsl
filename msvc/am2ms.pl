#!/usr/bin/perl

# This script is a mess.

open(LIBGSL, ">libgsl.dsp");
print LIBGSL &begin_project_lib("libgsl");
print LIBGSL &begin_target("libgsl");

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
    
    undef %AM;

    @lines = split("\n",$in);
    for (@lines) {
        s/#.*//g;
        ($var,$data) = split(/\s*=\s*/, $_, 2);
        $AM{$var} = $data;
#    print "$var IS $data\n";
    }

    @instlibs = parse_list($AM{lib_LTLIBRARIES});
    @noinstlibs = parse_list($AM{noinst_LTLIBRARIES});
    @gsllib = grep(/libgsl\.la/, @instlibs);
    @otherlibs = grep(!/libgsl\.la/, @instlibs);
    @ext_headers = parse_list($AM{pkginclude_HEADERS});
    @int_headers = parse_list($AM{noinst_HEADERS});
    @progs = parse_list($AM{bin_PROGRAMS}, $AM{check_PROGRAMS});

    push(@installed_headers, map("$dir/$_", @ext_headers));

    #print "== dir $dir ==\n" ;

    for $lib (@gsllib, @noinstlibs) {
        my $name = base($lib);
        #print "adding to libgsl ", $name, "\n";
        print LIBGSL &begin_group($name);
        print LIBGSL &add_files("libgsl",$dir, $name, 'inc', &list_sources($lib));
        print LIBGSL &add_files("libgsl",$dir, $name, 'exclude', grep(!/^test_/,@int_headers));
        print LIBGSL &add_files("libgsl",$dir, $name, 'exclude', @ext_headers);
        print LIBGSL &end_group();
    }

    for $lib (@otherlibs) {
        #print "write file for ", base($lib), "\n";
        my $name = base($lib);
        #print "adding to libgsl ", $name, "\n";
        open(LIB, ">$name.dsp");
        print LIB &begin_project_lib($name);
        print LIB &begin_target($name);
        #print LIB &begin_group($name);
        print LIB &add_files($name,$dir, $name, 'inc', &list_sources($lib));
        print LIB &add_files($name,$dir, $name, 'exclude', grep(!/^test_/,@int_headers));
        print LIB &add_files($name,$dir, $name, 'exclude', @ext_headers);
        #print LIB &end_group();
        print LIB &end_target();
        print LIB &end_project();
        close(LIB);
    }

    for $prog (@progs) {
        #print "put in tests $prog\n";
        @sources = &list_sources($prog);
        @ldadds = &list_ldadds($prog);
        my $name = $dir; 
        $name =~ s/.*\///; 
        $name =~ s/-/_/g;
        if ($prog =~ /test/) {
            $name = "${prog}_${name}";
            push (@tests, "$name");
        } else {
            $name = $prog;
        }
        open(TEST, ">$name.dsp");
        warn "$name $prog => ", @sources, "\n";
        print TEST &begin_project_app($name);
        print TEST &begin_target($name);
        #print TEST &begin_group("Source Files");
        print TEST &add_files("test", $dir, $name, 'inc', &list_sources($prog));
        print TEST &add_files("test", $dir, $name, 'exclude', grep(/^test_/,@int_headers));
        #print TEST &add_files("test", $dir, $name, 'exclude', @ext_headers);
        #print TEST &end_group("Source Files");
        print TEST &end_target();
        print TEST &end_project();
        close (TEST);
    }
}

print LIBGSL &end_target();
print LIBGSL &end_project();

#  Write GSL workspace

open(GSL, ">GSL.dsw");
print GSL &begin_workspace();
print GSL &add_workspace_project("libgsl", "libgsl.dsp");
print GSL &add_workspace_project("libgslcblas", "libgslcblas.dsp");
print GSL &end_workspace();
close GSL;

#  Write test workspace

open(TEST, ">GSLTESTS.dsw");
print TEST &begin_workspace();
print TEST &add_workspace_project("GSLTESTS", "GSLTESTS.dsp", @tests);
for $t (@tests) {
    print TEST &add_workspace_project($t, "$t.dsp");
}
print TEST &end_workspace();
close (TEST);

open(TEST, ">GSLTESTS.dsp");
print TEST &generic_project("GSLTESTS");
close(TEST);

# Write test batch files

open(TEST, ">MAKE_CHECK.bat");
print TEST &begin_check();
for $t (@tests) {
    print TEST &add_check($t);
}
print TEST &end_check();
close(TEST);

# Write script to copy headers into gsl directory

open(BATCH, ">COPY_GSL_HEADERS.bat");
for $h (@installed_headers) {
    $h =~ s#/#\\#g;  # convert to dos style path
    print BATCH "copy $h ..\\gsl\n" ;
}
close(BATCH);


######################################################################

sub parse_list {
    return split(' ', join(' ',@_));
}

sub list_sources {
    my ($f) = @_;
    $f =~ s/[\.\-]/_/g;
    return split(' ', $AM{"${f}_SOURCES"});
}

sub list_ldadds {
    my ($f) = @_;
    $f =~ s/[\.\-]/_/g;
    return split(' ', $AM{"${f}_LDADD"});
}

sub base {
    my ($f) = @_ ;
    $f =~ s/\.\w+$//;
    return $f;
}

######################################################################


sub begin_project_lib {
    my ($proj) = @_;

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
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MT /Za /W3 /GX /O2 /I "..\\msvc" /I "." /I ".." /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
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
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ  /c
# ADD CPP /nologo /MDd /Za /W3 /Gm /GX /ZI /Od  /I "..\\msvc" /I "." /I ".." /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ  /c
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

sub begin_project_app {
    my ($proj) = @_;

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
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MD /Za /W3 /GX /O2 /I "..\\msvc" /I "." /I ".." /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc0a /d "NDEBUG"
# ADD RSC /l 0xc0a /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib  kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 libgsl.lib libgslcblas.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386 /out:"bin\\$proj.exe" /libpath:"Release"

!ELSEIF  "\$(CFG)" == "$proj - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ  /c
# ADD CPP /nologo /MDd /Za /W3 /Gm /GX /ZI /Od /I "..\\msvc" /I "." /I ".." /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FD /GZ  /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0xc0a /d "_DEBUG"
# ADD RSC /l 0xc0a /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib  kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 libgsl.lib libgslcblas.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /out:"bin\\$proj.exe" /pdbtype:sept /libpath:"Debug"

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

sub add_files {
    my ($proj, $dir, $name, $flag, @files) = @_;
    my $file;
    my $y;
    $dir =~ s/\//\\/g; # convert to dos style path with \'s
    my $prop;
    $prop = "# PROP Exclude_From_Build 1" if $flag eq 'exclude';
    for $file (@files) {
    my $x = 
"# Begin Source File

SOURCE=$dir\\$file

$prop

!IF  \"\$(CFG)\" == \"$proj - Win32 Release\"

# PROP Intermediate_Dir \"Release\\$name\"

!ELSEIF  \"\$(CFG)\" == \"$proj - Win32 Debug\"

# PROP Intermediate_Dir \"Debug\\$name\"

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
    my ($name) = @_;
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
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""

!ELSEIF  "\$(CFG)" == "$name - Win32 Debug"

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
    my ($test) = @_;
    return <<"EOF";
\@echo running $test...
\@.\\bin\\$test.exe >> result.dat
EOF
}

sub end_check {
    return <<'EOF'
@echo Number of failures:
@FIND /V /C "PASS:" result.dat
@echo Full test log is in result.dat
EOF
}

######################################################################
