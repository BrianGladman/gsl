#!/usr/bin/perl

$/=undef;
$in = <>;

$in =~ s/\\\n/ /g;
$in =~ s/\n\n/\n/g;
$in =~ s/[ \t][ \t]+/ /g;
$in =~ s/^\s+//g;

@lines = split("\n",$in);
for (@lines) {
    next if /^\s*#/;
    ($var,$data) = split(/\s*=\s*/, $_, 2);
    $AM{$var} = $data;
#    print "$var IS $data\n";
}

# libraries are noinst_LTLIBRARIES lib_LTLIBRARIES
# installed headers are pkginclude_HEADERS
# headers are noinst_HEADERS
# includes are INCLUDES
# program bin_PROGRAMS check_PROGRAMS

@libs = parse_list($AM{noinst_LTLIBRARIES}, $AM{lib_LTLIBRARIES});
@ext_headers = parse_list($AM{pkginclude_HEADERS});
@int_headers = parse_list($AM{noinst_HEADERS});

for $lib (@libs) {
    print "# lib = $lib\n" ;
    @sources = &list_sources($lib);
    print "TEMPLATE = lib\n";
    print "CONFIG = warn_on release\n";
    print "SOURCES = ", join(" ", @sources), "\n" ; 
    print "INCLUDEPATH = ..\n";
    #print "INTERFACES = ", join(" ", @ext_headers), "\n" ; 
    print "HEADERS = ", join(" ", @ext_headers, @int_headers), "\n" ; 
    print "TARGET = ", base($lib), "\n" ;

}

@progs = parse_list($AM{bin_PROGRAMS}, $AM{check_PROGRAMS});

for $prog (@progs) {
    print "# app = $prog\n" ;
    @sources = &list_sources($prog);
    print "# SOURCES = ", join(" ", @sources), "\n" ; 
    @ldadds = &list_ldadds($prog);
    print "# needs libs ", join(", ", @ldadds), "\n" ; 
}


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
