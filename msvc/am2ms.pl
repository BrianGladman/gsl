#!/usr/bin/perl

print &begin_project("libgsl");
print &begin_target("libgsl");

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
        next if /^\s*#/;
        ($var,$data) = split(/\s*=\s*/, $_, 2);
        $AM{$var} = $data;
#    print "$var IS $data\n";
    }

    @instlibs = parse_list($AM{lib_LTLIBRARIES});
    @noinstlibs = parse_list($AM{noinst_LTLIBRARIES});
    @otherlibs = parse_list($AM{lib_LTLIBRARIES});
    @ext_headers = parse_list($AM{pkginclude_HEADERS});
    @int_headers = parse_list($AM{noinst_HEADERS});
    @progs = parse_list($AM{bin_PROGRAMS}, $AM{check_PROGRAMS});

    #print "== dir $dir ==\n" ;

    for $lib (@instlibs, @noinstlibs) {
        #print "adding to libgsl ", base($lib), "\n";
        print &begin_group(base($lib));
        print &add_files($dir, base($lib), 'inc', &list_sources($lib));
        print &add_files($dir, base($lib), 'exclude', grep(!/^test_/,@int_headers));
        print &add_files($dir, base($lib), 'exclude', @ext_headers);
        print &end_group();
    }

    
    for $lib (@otherlibs) {
        #print "write file for ", base($lib), "\n";
    }

    for $prog (@progs) {
        #print "put in tests $prog\n";
        @sources = &list_sources($prog);
        @ldadds = &list_ldadds($prog);
    }
}

print &end_target();
print &end_project();

# libraries are noinst_LTLIBRARIES lib_LTLIBRARIES
# installed headers are pkginclude_HEADERS
# headers are noinst_HEADERS
# includes are INCLUDES
# program bin_PROGRAMS check_PROGRAMS

# {
#     print "# lib = $lib\n" ;
#     print "TEMPLATE = lib\n";
#     print "SOURCES = ", join(" ", @sources), "\n" ; 
#     print "INCLUDEPATH = ..\n";
#     #print "INTERFACES = ", join(" ", @ext_headers), "\n" ; 
#     print "HEADERS = ", join(" ", @ext_headers, @int_headers), "\n" ; 
#     print "TARGET = ", base($lib), "\n" ;
# }
#     print "# app = $prog\n" ;
#     @sources = &list_sources($prog);
#     print "# SOURCES = ", join(" ", @sources), "\n" ; 
#     @ldadds = &list_ldadds($prog);
#     print "# needs libs ", join(", ", @ldadds), "\n" ; 
# }


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


sub begin_project {
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
# ADD CPP /nologo /W3 /GX /O2 /I "..\\msvc" /I "." /I ".." /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
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
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od  /I "..\\msvc" /I "." /I ".." /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ  /c
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
    my ($dir, $name, $flag, @files) = @_;
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

!IF  \"\$(CFG)\" == \"libgsl - Win32 Release\"

# PROP Intermediate_Dir \"Release\\$name\"

!ELSEIF  \"\$(CFG)\" == \"libgsl - Win32 Debug\"

# PROP Intermediate_Dir \"Debug\\$name\"

!ENDIF 

# End Source File
";
    $y .= $x;}
    return $y;
}

######################################################################
