#!/usr/bin/perl
use File::Basename;


$file = $ARGV[0];
$lib = basename($file); 
($name = $lib) =~ s/\..*$//;

print "LIBRARY $name\n";
print "DESCRIPTION \"GNU Scientific Library $lib\"\n";
print "EXPORTS\n";

open (LIB, "<$file") || die "can't open lib with nm: $!";
while (<LIB>) {
    chomp;
    my ($address, $type, $fn) = split(' ', $_);
    if ($type eq 'T') {
        print "\t$fn\n";
    } elsif ($type eq 'D') {
        print "\t$fn DATA\n";
    }
}
close(LIB);

