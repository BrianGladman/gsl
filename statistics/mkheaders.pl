#!/usr/bin/perl

# script to generate all the headers from gsl_statistics_int.h

open(MASTER,"<gsl_statistics_int.h") ;
while (<MASTER>) {
    s/^int /BASE / ;
    s/int data/BASE data/g ;
    s/int sorted_data/BASE sorted_data/g ;
    s/_int_/_BASE_/g ;
    s/_INT_/_UBASE_/g ;
    push(@lines,$_) ;
}
close(MASTER) ;

for $t ('double','float','long double', 
	'char', 'unsigned char', 
	'short', 'unsigned short', 
	'unsigned int', 
	'long', 'unsigned long')
{
    $l = $t ;
    $l =~ s/unsigned /u/ ;
    $l =~ s/ /_/;
    $u = "\U$l\E" ;
    $s = "$l\_" ; if ($t eq 'double') { $s = "" } ;

    open(FILE,">gsl_statistics_$l.h") ;

    @a = @lines ;
    for (@a) {
	s/_UBASE_/\_$u\_/g ;
	s/_BASE_/\_$s/g ;
	s/BASE/$t/g ;
    }
    print FILE @a ;
    close(FILE) ;
}
