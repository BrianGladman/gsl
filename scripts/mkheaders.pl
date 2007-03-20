#!/usr/bin/perl

# This script generates header files automatically
# from the "float" version of the header
# e.g mkheaders.pl gsl_statistics
# generates double, int, etc from gsl_statistics_float.h

$dir = $ARGV[0];
$f = "${dir}_float.h" ;
$x = "float" ;
#$suff = ".new";

die "can't find $f" if ! -e $f;

open(FILE,"<$f") ; @lines = <FILE> ; close(FILE) ;

&substitute("double") ;
&substitute("long double") ;
&substitute("int") ;
&substitute("char") ;
&substitute("short") ;
&substitute("long") ;
&substitute("unsigned int") ;
&substitute("unsigned char") ;
&substitute("unsigned short") ;
&substitute("unsigned long") ;

$f = "${dir}_complex_float.h" ;
$x = "float" ;
exit if ! -e $f ;
open(FILE,"<$f") ; @lines = <FILE> ; close(FILE) ;
&substitute("double") ;
&substitute("long double") ;

sub substitute {
    my ($t) = @_ ;
    
    $s = $t ;
    if ($t eq "double") {
        $s = "" ;
        $u = $t;
    } else {
        $s =~ s/^unsigned /u/;
        $s =~ s/ /_/g;
        $u = $s;
        $s = "_$s" ;
    }

    ($new = $f) =~ s/$x/$u/;
    die if $new eq $f ;
    @temp = @lines ;
    print $new , "\n";
    open (FILE,">$new$suff");
    for (@temp) {
        $skip = 1 if /DEFINED FOR FLOATING POINT TYPES ONLY/i && $t !~ /float|double/;
        $skip = 0, next if /END OF FLOATING POINT TYPES/i && $skip;
        next if $skip ;

        s/_$x\.h/_$u.h/g;
        s/_$x/$s/g;
        s/$x/$t/g;
        s/\U${x}_H/\U${u}_H/g;
        s/_\U$x/\U$s/g;
        print FILE $_ ;
    }
    close(FILE);
}
