#!/usr/bin/perl
use File::Basename;


$file = $ARGV[0];
$lib = basename($file); 
$lib =~ s/^lib//;
($name = $lib) =~ s/\..*$//;

print "LIBRARY $name\n";
print "DESCRIPTION \"GNU Scientific Library ${name}.lib\"\n";
print "EXPORTS\n";

open (LIB, "<$file") || die "can't open file: $!";
while (<LIB>) {
    chomp;
    my ($address, $type, $fn) = split(' ', $_);
    if ($type eq 'T') {
        print "\t$fn\n";
    } elsif ($type eq 'D' or $type eq 'R' or $type eq 'C') {
        print "\t$fn DATA\n";
    }
}
close(LIB);

